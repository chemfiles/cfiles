// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <fstream>

#include <docopt/docopt.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include "Rotcf.hpp"
#include "Autocorrelation.hpp"
#include "warnings.hpp"

using namespace chemfiles;

static const char OPTIONS[] =
R"(Compute rotation correlation dynamic for arbitrary bonds and molecules. The
bonds and molecules to use are specified using chemfiles selection language.
This analysis does not support changes in the topology or the matched atoms
during the simulation.

For more information about chemfiles selection language, please see
http://chemfiles.org/chemfiles/latest/selections.html

Usage:
  cfiles rotcf [options] <trajectory>
  cfiles rotcf (-h | --help)

Examples:
  cfiles rotcf water.xyz --cell 15:15:25
  cfiles rotcf input.pdb -s "bonds: type(#1) O and type(#2) H"

Options:
  -h --help                     show this help
  -o <file>, --output=<file>    write result to <file>. This default to the
                                trajectory file name with the `.rotcf.dat`
                                extension.
  --format=<format>             force the input file format to be <format>
  -t <path>, --topology=<path>  alternative topology file for the input
  --topology-format=<format>    use <format> as format for the topology file
  --guess-bonds                 guess the bonds in the input
  -c <cell>, --cell=<cell>      alternative unit cell. <cell> format is one of
                                <a:b:c:α:β:γ> or <a:b:c> or <a>. 'a', 'b' and
                                'c' are in angstroms, 'α', 'β', and 'γ' are in
                                degrees.
  --steps=<steps>               steps to use from the input. <steps> format
                                is <start>:<end>[:<stride>] with <start>, <end>
                                and <stride> optional. The used steps goes from
                                <start> to <end> (excluded) by steps of
                                <stride>. The default values are 0 for <start>,
                                the number of steps for <end> and 1 for
                                <stride>.
  --selection=<sel>, -s <sel>   selection to use for the donors. This must be a
                                selection of size 2 [default: bonds: all]
)";

static Rotcf::Options parse_options(int argc, const char* argv[]) {
    auto options_str = command_header("rotcf", Rotcf().description()) + "\n";
    options_str += "Guillaume Fraux <guillaume@fraux.fr>\n";
    options_str += OPTIONS;
    auto args = docopt::docopt(options_str, {argv, argv + argc}, true, "");

    Rotcf::Options options;
    options.trajectory = args.at("<trajectory>").asString();
    options.guess_bonds = args.at("--guess-bonds").asBool();

    options.selection = args.at("--selection").asString();

    if (args.at("--output")) {
        options.outfile = args.at("--output").asString();
    } else {
        options.outfile = options.trajectory + ".rotcf.dat";
    }

    if (args.at("--steps")) {
        options.steps = steps_range::parse(args.at("--steps").asString());
    }

    if (args.at("--format")) {
        options.format = args.at("--format").asString();
    }

    if (args.at("--topology")) {
        if (options.guess_bonds) {
            throw CFilesError("Can not use both '--topology' and '--guess-bonds'");
        }
        options.topology = args.at("--topology").asString();
    }

    if (args.at("--topology-format")) {
        if (options.topology == "") {
            throw CFilesError("Can not use '--topology-format' without a '--topology'");
        }
        options.topology_format = args.at("--topology-format").asString();
    }

    if (args.at("--cell")) {
        options.custom_cell = true;
        options.cell = parse_cell(args.at("--cell").asString());
    }

    return options;
}

std::string Rotcf::description() const {
    return "rotation correlation dynamic for arbitrary bonds and molecules";
}

int Rotcf::run(int argc, const char* argv[]) {
    auto options = parse_options(argc, argv);

    auto selection = Selection(options.selection);
    if (selection.size() != 2) {
        throw CFilesError("Selection must have a size of 2 (either bonds: or pairs:)");
    }

    auto trajectory = Trajectory(options.trajectory, 'r', options.format);
    if (options.custom_cell) {
        trajectory.set_cell(options.cell);
    }

    if (options.topology != "") {
        trajectory.set_topology(options.topology, options.topology_format);
    }

    auto frame = trajectory.read_step(options.steps.first());
    if (options.guess_bonds) {
        frame.guess_bonds();
    }

    auto matched = selection.evaluate(frame);
    if (matched.empty()) {
        warn("no matching atom in the first frame");
        return 0;
    }

    auto vectors = std::vector<std::vector<Vector3D>>(matched.size());
    for (auto step: options.steps) {
        if (step >= trajectory.nsteps()) {
            break;
        }
        auto frame = trajectory.read_step(step);
        auto positions = frame.positions();
        for (size_t i=0; i<matched.size(); i++) {
            auto& match = matched[i];
            assert(match.size() == 2);

            auto rij = frame.cell().wrap(positions[match[0]] - positions[match[1]]);
            rij /= rij.norm();
            vectors[i].push_back(rij);
        }
    }

    // Following GROMACS, we compute the P2 autocorrelation using 6 different
    // FFT:
    //
    // C2(t) = <P2(u(0) ⋅ u(t))>
    //       = <1/2 (3 * (u(0) ⋅ u(t))^2 - 1)>
    //       = <1/2 (3 * cos^2(θ) - 1)>
    //       = 3/2 (<x^2> + <y^2> + <z^2> + 2<xy> + 2<xz> + 2<yz>) - 1/2

    // Accessing vectors[0] is fine, as we already exited if no atoms matched
    // the selection.
    auto used_steps = vectors[0].size();
    auto result = std::vector<float>(used_steps / 2, 0.0);

    auto do_correlation = [&](size_t i, size_t j) {
        auto correlator = Autocorrelation(used_steps);
        for (auto& vector: vectors) {
            auto squares = std::vector<float>(vector.size());
            for (size_t step=0; step<used_steps; step++) {
                squares[step] = vector[step][i] * vector[step][j];
            }
            correlator.add_timeserie(std::move(squares));
        }
        correlator.normalize();
        auto& correlation = correlator.get_result();
        auto factor = i == j ? 1.5 : 3.0;
        for (size_t step=0; step<result.size(); step++) {
            result[step] += factor * correlation[step];
        }
    };

    do_correlation(0, 0);
    do_correlation(1, 1);
    do_correlation(2, 2);
    do_correlation(0, 1);
    do_correlation(0, 2);
    do_correlation(1, 2);

    for (size_t i=0; i<result.size(); i++) {
        result[i] -= 0.5;
    }

    std::ofstream output(options.outfile, std::ios::out);
    if (!output.is_open()) {
        throw CFilesError("Could not open the '" + options.outfile + "' file.");
    }
    fmt::print(output, "# rotation correlation for \"{}\" in {}\n", options.selection, options.trajectory);
    fmt::print(output, "# step value\n");

    for (size_t i=0; i<result.size(); i++) {
        fmt::print(output, "{} {}\n", i * options.steps.stride(), result[i]);
    }

    return 0;
}
