// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#include <docopt/docopt.h>
#include <sstream>
#include <fstream>

#include "HBonds.hpp"
#include "Errors.hpp"
#include "utils.hpp"
#include "geometry.hpp"

using namespace chemfiles;

static const double pi = 3.141592653589793238463;
static const char OPTIONS[] =
R"(Compute list of hydrogen bonds along a trajectory. A selection can be specified using the 
chemfiles selection language. It is possible to provide an alternative unit cell or topology 
for the trajectory file if they are not defined in the trajectory format. Hydrogen bonds criteria
can be specified (donor-acceptor distance and donor-acceptor-H angle)

For more information about chemfiles selection language, please see
http://chemfiles.github.io/chemfiles/latest/selections.html

Usage:
  cfiles hbonds [options] <trajectory>
  cfiles hbonds (-h | --help)

Examples:
  cfiles hbonds --cell=28 --guess-bonds water.xyz water.pdb
  cfiles hbonds butane.pdb butane.nc --wrap
  cfiles hbonds methane.xyz --cell 15:15:25 --guess-bonds --points=150
  cfiles hbonds result.xtc --topology=initial.mol --topology-format=PDB out.nc
  cfiles hbonds in.zeo out.mol --input-format=XYZ --output-format=PDB

Options:
  -h --help                     show this help
  -o <file>, --output=<file>    write result to <file>. This default to the 
                                trajectory file name with the `.hb` extension.
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
                                and <stride> optional. Default is to use all
                                steps from the input; starting at 0, ending at
                                the last step, and with a stride of 1.
  -s <sel>, --selection=<sel>   selection to use for the atoms. This must be a
                                selection of size 3 (for angles) or 4 (for
                                dihedral angles) [default: angles: all]
  -p <par>, --parameters=<par>  parameters to use for the hydrogen bond. <par>
                                format is <d:α> where 'd' is the donor-acceptor 
                                maximum distance in angstroms and 'α' is the
                                donor-acceptor-hydrogen maximum angle in degrees.
                                [default: 3.0:30.0]
)";

static HBonds::Options parse_options(int argc, const char* argv[]) {
    auto options_str = command_header("hbonds", HBonds().description()) + "\n";
    options_str += OPTIONS;
    auto args = docopt::docopt(options_str, {argv, argv + argc}, true, "");

    HBonds::Options options;
    options.trajectory = args["<trajectory>"].asString();
    options.guess_bonds = args.at("--guess-bonds").asBool();
    options.selection = args["--selection"].asString();

    if (args["--output"]){
        options.outfile = args["--output"].asString();
    } else {
        options.outfile = options.trajectory + ".hb";
    }

    if (args.at("--steps")) {
        options.steps = steps_range::parse(args.at("--steps").asString());
    }

    if (args["--format"]){
        options.format = args["--format"].asString();
    }

    if (args["--topology"]){
        if (options.guess_bonds) {
            throw CFilesError("Can not use both '--topology' and '--guess-bonds'");
        }
        options.topology = args["--topology"].asString();
    }

    if (args["--topology-format"]){
        if (options.topology == "") {
            throw CFilesError("Useless '--topology-format' without '--topology'");
        }
        options.topology_format = args["--topology-format"].asString();
    }

    if (args["--cell"]) {
        options.custom_cell = true;
	options.cell = parse_cell(args["--cell"].asString());
	}

    if (args["--parameters"]) {
	auto splitted = split(args["--parameters"].asString(), ':');
	if (splitted.size() == 2) {
	    options.parameters[0] = stod(splitted[0]);
	    options.parameters[1] = stod(splitted[1]);
	} else {
            throw CFilesError(
                "custom parameters should be specified as 'd:α'"
            );
    	}
    }	

    return options;
}


std::string HBonds::description() const {
    return "compute hydrogen bonds network";
}

int HBonds::run(int argc, const char* argv[]) {
    auto options = parse_options(argc, argv);

    auto infile = Trajectory(options.trajectory, 'r', options.format);
    std::ofstream outfile(options.outfile, std::ios::out);
    if (outfile.is_open()) {
        outfile << "#Hydrogen bond network in trajectory " << options.trajectory << std::endl;
        outfile << "# Selection: " << options.selection << std::endl;
        outfile << "# Criteria:" << std::endl;
        outfile << "# donor-acceptor distance < " << options.parameters[0] << " angstroms" << std::endl;
        outfile << "# donor-acceptor-H angle < " << options.parameters[1] << " degrees" << std::endl;
    } else {
        throw CFilesError("Could not open the '" + options.outfile + "' file.");
    }

    if (options.custom_cell) {
        infile.set_cell(options.cell);
    }

    if (options.topology != "") {
        infile.set_topology(options.topology, options.topology_format);
    }

    for (auto step: options.steps) {
        if (step >= infile.nsteps()) {
            break;
        }
        auto frame = infile.read_step(step);
        if (options.guess_bonds) {
            frame.guess_topology();
        }

        outfile << "# Frame: " << step << std::endl;        

        auto positions = frame.positions();
        auto cell = frame.cell();

        auto matched = selectionAcceptor_.evaluate(frame);
        outfile << matched.size() << std::endl;
        for (auto match: matched) {
            outfile << "entre dans match" << std::endl;
            assert(match.size() == 2);

            auto acceptor = match[0];
            auto hydrogen = match[1];
            auto donors = selectionDonor_.evaluate(frame);
            for (auto donor: donors) {
                auto rad = cell.wrap(positions[donor[0]] - positions[acceptor]);
                auto distance = norm(rad); 
                auto rah = cell.wrap(positions[hydrogen] - positions[acceptor]);
                auto theta = angle(rad, rah)*pi/180;
                outfile << distance << "    " << theta << std::endl;
                if (distance < options.parameters[0] && theta < options.parameters[1]) {
                    outfile << "LH " << distance << "    " << theta << std::endl;
                }
            }
        }
    }

    return 0;
}
