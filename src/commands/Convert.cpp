// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <docopt/docopt.h>
#include <sstream>
#include <set>

#include "Convert.hpp"
#include "Errors.hpp"
#include "utils.hpp"

using namespace chemfiles;

static const char OPTIONS[] =
R"(Convert trajectories from one format to another, automatically guessing the
format to use based on the files extension. It is possible to force a specific
input or output file format, and to specify an alternative unit cell or topology
for the input file if they are not defined in the input format.
One may write only a part of the input file by defining a selection using
the chemfiles selection language.

For more information about chemfiles selection language, please see
http://chemfiles.org/chemfiles/latest/selections.html

Usage:
  cfiles convert [options] <input> <output>
  cfiles convert (-h | --help)

Examples:
  cfiles convert --cell=28 --guess-bonds water.xyz water.pdb
  cfiles convert butane.pdb butane.nc --wrap
  cfiles convert methane.xyz --cell 15:15:25 --guess-bonds --points=150
  cfiles convert result.xtc --topology=initial.mol --topology-format=PDB out.nc
  cfiles convert in.zeo out.mol --input-format=XYZ --output-format=PDB

Options:
  -h --help                     show this help
  --input-format=<format>       force the input file format to be <format>
  --output-format=<format>      force the output file format to be <format>
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
  --wrap                        rewrap the particles inside the unit cell
  -s <sel>, --selection=<sel>   selection to use for the output file
                                [default: atoms: all]
)";

static Convert::Options parse_options(int argc, const char* argv[]) {
    auto options_str = command_header("convert", Convert().description());
    options_str += "Guillaume Fraux <guillaume@fraux.fr>\n\n";
    options_str += OPTIONS;
    auto args = docopt::docopt(options_str, {argv, argv + argc}, true, "");

    Convert::Options options;
    options.infile = args.at("<input>").asString();
    options.outfile = args.at("<output>").asString();
    options.guess_bonds = args.at("--guess-bonds").asBool();
    options.wrap = args.at("--wrap").asBool();
    options.selection = args.at("--selection").asString();

    if (args.at("--steps")) {
        options.steps = steps_range::parse(args.at("--steps").asString());
    }

    if (args.at("--input-format")){
        options.input_format = args.at("--input-format").asString();
    }

    if (args.at("--output-format")){
        options.output_format = args.at("--output-format").asString();
    }

    if (args.at("--topology")){
        if (options.guess_bonds) {
            throw CFilesError("Can not use both '--topology' and '--guess-bonds'");
        }
        options.topology = args.at("--topology").asString();
    }

    if (args.at("--topology-format")){
        if (options.topology == "") {
            throw CFilesError("Can not use '--topology-format' without a '--topology'");
        }
        options.topology_format = args["--topology-format"].asString();
    }

    if (args.at("--cell")) {
        options.custom_cell = true;
        options.cell = parse_cell(args.at("--cell").asString());
    }

    return options;
}


std::string Convert::description() const {
    return "convert trajectories between formats";
}

int Convert::run(int argc, const char* argv[]) {
    auto options = parse_options(argc, argv);

    auto infile = Trajectory(options.infile, 'r', options.input_format);
    auto outfile = Trajectory(options.outfile, 'w', options.output_format);

    if (options.custom_cell) {
        infile.set_cell(options.cell);
    }

    if (options.topology != "") {
        infile.set_topology(options.topology, options.topology_format);
    }

    selection_ = Selection(options.selection);
    for (auto step: options.steps) {
        if (step >= infile.nsteps()) {
            break;
        }
        auto frame = infile.read_step(step);

        if (options.guess_bonds) {
            frame.guess_topology();
        }

        if (options.wrap) {
            auto positions = frame.positions();
            auto cell = frame.cell();

            for (auto& position: positions) {
                position = cell.wrap(position);
            }
        }

        if (options.selection != "") {
            auto matched = selection_.evaluate(frame);

            std::set<size_t> keep;
            for (auto match: matched) {
                for (size_t i = 0; i < match.size(); i++) {
                    keep.insert(match[i]);
                }
            }

            std::vector<size_t> remove;
            remove.reserve(frame.natoms() - keep.size());
            for (size_t i = 0; i < frame.natoms(); ++i) {
                auto search = keep.find(i);
                if (search == keep.end()) { // element not found
                    remove.push_back(i);
                }
            }

            // deleting an atom in the frame shifts all the indexes after it
            // so we need to iterate in reverse order to delete the good atoms
            for (auto i: reverse(remove)) {
                frame.remove(i);
            }
        }

        outfile.write(frame);
    }

    return 0;
}
