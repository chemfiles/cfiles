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
R"(Compute list of hydrogen bonds along a trajectory. Selections for the acceptor and donor atoms 
can be specified using the chemfiles selection language. It is possible to provide an alternative 
unit cell or topology for the trajectory file if they are not defined in the trajectory format. 
Hydrogen bonds are defined as electrostatic attraction between two polar groups: the acceptor group 
is a hydrogen atom covalently bound to an electronegative atom (usually O, N, F) while 
the donor group is another highly electronegative atom. 
The criteria used depend on a maximum donor-acceptor distance and a maximum donor-acceptor-H angle.
Hydrogen bonds criteria can be specified.

For more information about chemfiles selection language, please see
http://chemfiles.org/chemfiles/latest/selections.html

Usage:
  cfiles hbonds [options] <trajectory>
  cfiles hbonds (-h | --help)

Examples:
  cfiles hbonds water.xyz --cell 15:15:25 --guess-bonds 
  cfiles hbonds in.pdb --acceptors=="bonds: type(#1) == O and type(#2) == H"
  cfiles hbonds protein.pdb --donors=="atoms: type N" -p 2.5:20.0

Options:
  -h --help                     show this help
  -o <file>, --output=<file>    write result to <file>. This default to the 
                                trajectory file name with the `_hb.dat` extension.
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
  --acceptors=<sel>             selection to use for the acceptors. This must be a
                                selection of size 2 and type bonds The second atom
				should be the hydrogen atom.
                                [default: bonds: type(#2) == H]
  --donors=<sel>                selection to use for the donors. This must be a
                                selection of size 1.
                                [default: atoms: type O or type N or type F]
  --distance=<dist>             distance criterion to use for the hydrogen bond detection.
                                'dist' is the donor-acceptor maximum distance in angstroms.
                                [default: 3.0]
  --angle=<angle>               angle criterion to use for the hydrogen bond detection.
                                'angle' is the donor-acceptor-hydrogen maximum angle in degrees.
                                [default: 30.0]
)";

static HBonds::Options parse_options(int argc, const char* argv[]) {
    auto options_str = command_header("hbonds", HBonds().description()) + "\n";
    options_str += "Laura Scalfi <laura.scalfi@ens.fr>\n";
    options_str += OPTIONS;
    auto args = docopt::docopt(options_str, {argv, argv + argc}, true, "");

    HBonds::Options options;
    options.trajectory = args.at("<trajectory>").asString();
    options.guess_bonds = args.at("--guess-bonds").asBool();

    options.selectionAcceptor = args.at("--acceptors").asString();
    options.selectionDonor = args.at("--donors").asString();

    if (args.at("--output")){
        options.outfile = args.at("--output").asString();
    } else {
        options.outfile = options.trajectory + "_hb.dat";
    }

    if (args.at("--steps")) {
        options.steps = steps_range::parse(args.at("--steps").asString());
    }

    if (args.at("--format")){
        options.format = args.at("--format").asString();
    }

    if (args.at("--topology")){
        if (options.guess_bonds) {
            throw CFilesError("Can not use both '--topology' and '--guess-bonds'");
        }
        options.topology = args.at("--topology").asString();
    }

    if (args.at("--topology-format")){
        if (options.topology == "") {
            throw CFilesError("Useless '--topology-format' without '--topology'");
        }
        options.topology_format = args.at("--topology-format").asString();
    }

    if (args.at("--cell")) {
        options.custom_cell = true;
        options.cell = parse_cell(args.at("--cell").asString());
	}

    options.distance = 3.0;
    options.angle = 30.0;
    if (args.at("--distance")) {
        try {
            options.distance = std::stod(args.at("--distance").asString());
        } catch(std::invalid_argument) {
            throw CFilesError("Invalid argument for distance criterion");
        }
    }	
    if (args.at("--angle")) {
        try {
            options.angle = std::stod(args.at("--angle").asString())*pi/180;
        } catch(std::invalid_argument) {
            throw CFilesError("Invalid argument for angle criterion");
        }
    }	

    return options;
}


std::string HBonds::description() const {
    return "compute hydrogen bonds network";
}

int HBonds::run(int argc, const char* argv[]) {
    auto options = parse_options(argc, argv);

    selectionAcceptor_ = Selection(options.selectionAcceptor);
    if (selectionAcceptor_.size() != 2) {
        throw CFilesError("Can not use a selection for acceptors with size different than 2.");
    }

    selectionDonor_ = Selection(options.selectionDonor);
    if (selectionDonor_.size() != 1) {
        throw CFilesError("Can not use a selection for donors with size larger than 1.");
    }

    auto infile = Trajectory(options.trajectory, 'r', options.format);
    std::ofstream outfile(options.outfile, std::ios::out);
    if (outfile.is_open()) {
        outfile << "#Hydrogen bond network in trajectory " << options.trajectory << std::endl;
        outfile << "# Selection: acceptors: " << options.selectionAcceptor << " and donors: " << options.selectionDonor << std::endl;
        outfile << "# Criteria:" << std::endl;
        outfile << "# donor-acceptor distance < " << options.distance << " angstroms" << std::endl;
        outfile << "# donor-acceptor-H angle < " << options.angle*180/pi << " degrees" << std::endl;
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
        for (auto match: matched) {
            assert(match.size() == 2);

            auto acceptor = match[0];
            auto hydrogen = match[1];
            auto donors = selectionDonor_.list(frame);
            for (auto donor: donors) {
                if (donor != acceptor && frame.topology()[donor].type() != "H") {
                    auto r_ad = cell.wrap(positions[donor] - positions[acceptor]);
                    auto distance = norm(r_ad); 
                    auto r_ah = cell.wrap(positions[hydrogen] - positions[acceptor]);
                    auto theta = angle(r_ad, r_ah);
                    if (distance < options.distance && theta < options.angle) {
                        outfile << frame.topology()[donor].type() << donor << "   " << frame.topology()[acceptor].type() << acceptor << "   " << frame.topology()[hydrogen].type() << hydrogen << "  : " << distance << "    " << theta*180/pi << std::endl;
                    }
                }
            }
        }
    }

    return 0;
}
