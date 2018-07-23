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
  --wrap                        rewrap the particles matching the wrapping
                                selection inside the unit cell
  --wrap-selection=<self>       selection of atoms to wrap inside the cell
                                [default: all]
  --center                      translate all the atoms to place the center of
                                mass of the corresponding selection at the
                                origin. If both --wrap and --center are used,
                                the particles are wrapped first, and then
                                centered
  --center-selection=<self>     selection of atoms to use to compute the center
                                of mass to center inside the cell [default: all]
  -s <sel>, --selection=<sel>   selection to use for the output file
                                [default: all]
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
    options.wrap_selection = args.at("--wrap-selection").asString();
    options.center = args.at("--center").asBool();
    options.center_selection = args.at("--center-selection").asString();
    options.selection = args.at("--selection").asString();

    if (options.wrap_selection != "all" && !options.wrap) {
        throw cfiles_error("'--wrap-selection' without --wrap does nothing");
    }

    if (options.center_selection != "all" && !options.center) {
        throw cfiles_error("'--center-selection' without --center does nothing");
    }

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
            throw cfiles_error("Can not use both '--topology' and '--guess-bonds'");
        }
        options.topology = args.at("--topology").asString();
    }

    if (args.at("--topology-format")){
        if (options.topology == "") {
            throw cfiles_error("Can not use '--topology-format' without a '--topology'");
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

    auto selection = Selection(options.selection);
    auto wrap_sel = Selection(options.wrap_selection);
    if (wrap_sel.size() != 1) {
        throw cfiles_error("the wrapping selection should act on atoms");
    }
    auto center_sel = Selection(options.center_selection);
    if (center_sel.size() != 1) {
        throw cfiles_error("the center selection should act on atoms");
    }
    for (auto step: options.steps) {
        if (step >= infile.nsteps()) {
            break;
        }
        auto frame = infile.read_step(step);

        if (options.guess_bonds) {
            frame.guess_bonds();
        }

        if (options.wrap) {
            auto positions = frame.positions();
            auto cell = frame.cell();

            if (options.wrap_selection != "all") {
                for (auto i: wrap_sel.list(frame)) {
                    positions[i] = cell.wrap(positions[i]);
                }
            } else {
                for (auto& position: positions) {
                    position = cell.wrap(position);
                }
            }
        }

        if (options.center) {
            auto positions = frame.positions();
            double total_mass = 0.0;
            auto com = Vector3D();
            if (options.center_selection != "all") {
                for (size_t i: center_sel.list(frame)) {
                    auto mass = frame[i].mass();
                    com = com + mass * positions[i];
                    total_mass += mass;
                }
            } else {
                for (size_t i=0; i<frame.size(); i++) {
                    auto mass = frame[i].mass();
                    com = com + mass * positions[i];
                    total_mass += mass;
                }
            }
            com = com / total_mass;

            for (auto& position: positions) {
                position = position - com;
            }
        }

        if (options.selection != "all") {
            auto matched = selection.evaluate(frame);

            std::set<size_t> keep;
            for (auto match: matched) {
                for (size_t i = 0; i < match.size(); i++) {
                    keep.insert(match[i]);
                }
            }

            std::vector<size_t> remove;
            remove.reserve(frame.size() - keep.size());
            for (size_t i = 0; i < frame.size(); ++i) {
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
