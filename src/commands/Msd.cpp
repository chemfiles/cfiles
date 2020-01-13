// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <docopt/docopt.h>
#include <fstream>
#include <numeric>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include "Msd.hpp"
#include "Autocorrelation.hpp"
#include "Errors.hpp"
#include "utils.hpp"
#include "warnings.hpp"

using namespace chemfiles;

static const char OPTIONS[] =
R"(Compute mean square distance of an atom or a group of atom since the first
step of a trajectory. The resulting graph should be linear after a while, and
can be used to extract diffusion coefficient D for movement in d dimensions:
    <[r(t) - r(0)]^2> = 2 * d * D * t

Usage:
  cfiles msd [options] <trajectory>
  cfiles msd (-h | --help)

Examples:
  cfiles msd file.pdb -o msd.dat
  cfiles msd water.xyz --cell 15:15:25 --unwrap
  cfiles msd trajectory.nc --topology topol.pdb --selection "name Li"

Options:
  -h --help                     show this help
  -o <file>, --output=<file>    write result to <file>. This default to the
                                trajectory file name with the `.msd.dat`
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
  --selection=<sel>             selection of atoms to use when computing the
                                mean square distance. The selection should
                                always return the same atoms in the same order.
                                [default: all]
  --unwrap                      undo periodic boundary condition wrapping,
                                placing atoms back outside of the box
)";

static MSD::Options parse_options(int argc, const char* argv[]) {
    auto options_str = command_header("msd", MSD().description()) + "\n";
    options_str += "Guillaume Fraux <guillaume@fraux.fr>\n\n";
    options_str += OPTIONS;
    auto args = docopt::docopt(options_str, {argv, argv + argc}, true, "");

    MSD::Options options;
    options.trajectory = args.at("<trajectory>").asString();
    options.guess_bonds = args.at("--guess-bonds").asBool();

    options.selection = args.at("--selection").asString();

    if (args.at("--output")) {
        options.outfile = args.at("--output").asString();
    } else {
        options.outfile = options.trajectory + ".msd.dat";
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

    options.unwrap = args.at("--unwrap").asBool();

    return options;
}

std::string MSD::description() const {
    return "compute average mean square distance for a group of atoms";
}

int MSD::run(int argc, const char* argv[]) {
    auto options = parse_options(argc, argv);

    auto selection = Selection(options.selection);
    if (selection.size() != 1) {
        throw CFilesError("Can not use a selection with size larger than 1.");
    }

    std::ofstream outfile(options.outfile, std::ios::out);
    if (!outfile.is_open()) {
        throw CFilesError("Could not open the '" + options.outfile + "' file.");
    }
    fmt::print(outfile, "# Mean Square Deviation in {}\n", options.trajectory);
    fmt::print(outfile, "# For atoms '{}'\n", options.selection);

    auto trajectory = Trajectory(options.trajectory, 'r', options.format);
    if (options.custom_cell) {
        trajectory.set_cell(options.cell);
    }

    if (options.topology != "") {
        trajectory.set_topology(options.topology, options.topology_format);
    }

    // Pre-allocate memory to store the positions of each atom at each time step
    auto frame = trajectory.read_step(options.steps.first());
    if (options.guess_bonds) {
        frame.guess_bonds();
    }
    auto natoms = selection.list(frame).size();
    auto nsteps = options.steps.count(trajectory.nsteps());

    auto positions = std::vector<std::array<std::vector<float>, 3>>(natoms);
    for (size_t atom=0; atom<natoms; atom++) {
        positions[atom][0] = std::vector<float>(nsteps, 0.0);
        positions[atom][1] = std::vector<float>(nsteps, 0.0);
        positions[atom][2] = std::vector<float>(nsteps, 0.0);
    }

    // First, extract all the positions we need
    size_t current_step = 0;
    auto previous_frame = std::move(frame);
    for (auto step: options.steps) {
        if (step >= trajectory.nsteps()) {
            break;
        }
        auto frame = trajectory.read_step(step);
        if (options.guess_bonds) {
            frame.guess_bonds();
        }

        auto matched = selection.list(frame);
        if (matched.size() != natoms) {
            throw CFilesError(fmt::format(
                "the number of atoms matched by '{}' changed from {} to {} since the first step",
                options.selection, natoms, matched.size()
            ));
        }

        auto current_positions = frame.positions();
        auto previous_positions = previous_frame.positions();

        auto cell = frame.cell().matrix();
        auto prev_cell = previous_frame.cell().matrix();
        auto cell_inv = cell;
        auto prev_cell_inv = prev_cell;

        if (options.unwrap) {
            if (frame.cell().shape() == UnitCell::INFINITE) {
                throw CFilesError("can not unwrap in infinite unit cell");
            }
            cell_inv = cell.invert();
            prev_cell_inv = prev_cell.invert();
        }

        for (size_t atom=0; atom<natoms; atom++) {
            auto& current = current_positions[matched[atom]];

            if (options.unwrap) {
                auto curr_frac = cell_inv * current;
                auto prev_frac = prev_cell_inv * previous_positions[matched[atom]];
                auto delta = curr_frac - prev_frac;

                delta[0] -= round(delta[0]);
                delta[1] -= round(delta[1]);
                delta[2] -= round(delta[2]);

                current = cell * (prev_frac + delta);
            }

            positions[atom][0][current_step] = current[0];
            positions[atom][1][current_step] = current[1];
            positions[atom][2][current_step] = current[2];
        }

        current_step++;
        previous_frame = std::move(frame);
    }

    // We want to compute <[r(t) - r(0)]^2> where <...> denotes average on the
    // time origins and on the atoms. To do so, we separate the above expression
    // into <r(t)^2 + r(0)^2> - 2 <r(t) * r(0)>. The two first terms can be
    // computed directly, and the last one through the autocorrelation
    // framework.
    auto msd = std::vector<double>(nsteps, 0.0);

    // Start with the <r(t)^2 + r(0)^2> term
    for (size_t atom=0; atom<natoms; atom++) {
        auto rsq = std::vector<double>(nsteps, 0.0);
        for (size_t step=0; step<nsteps; step++) {
            auto xx = positions[atom][0][step] * positions[atom][0][step];
            auto yy = positions[atom][1][step] * positions[atom][1][step];
            auto zz = positions[atom][2][step] * positions[atom][2][step];

            rsq[step] = xx + yy + zz;
        }

        auto sum_rsq = 2 * std::accumulate(rsq.begin(), rsq.end(), 0.0);
        msd[0] += sum_rsq;

        double cum_sum = 0;
        double cum_sum_reverse = 0;
        for (size_t step=1; step<nsteps; step++) {
            cum_sum += rsq[step - 1];
            cum_sum_reverse += rsq[nsteps - step];
            msd[step] += (sum_rsq - cum_sum - cum_sum_reverse) / (nsteps - step);
        }
    }

    for (size_t step=0; step<nsteps; step++) {
        msd[step] /= natoms;
    }

    // compute the autocorrelation part
    auto correlation = Autocorrelation(nsteps);
    for (size_t atom=0; atom<natoms; atom++) {
        correlation.add_timeserie(std::move(positions[atom][0]));
        correlation.add_timeserie(std::move(positions[atom][1]));
        correlation.add_timeserie(std::move(positions[atom][2]));
    }
    correlation.normalize();

    auto& correlated = correlation.get_result();
    for (size_t step=1; step<nsteps; step++) {
        // the factor 3 is here because the correlation was normalized by
        // 3 * natoms (the total number of time series it got), but we need it
        // normalized by natoms only.
        msd[step] += -2 * 3 * correlated[step];
    }

    for (size_t step=1; step<nsteps / 2; step++) {
        fmt::print(outfile, "{} {}\n", step * options.steps.stride(), msd[step]);
    }

    return 0;
}
