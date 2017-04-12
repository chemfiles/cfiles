// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

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
Hydrogen bonds are defined as electrostatic attraction between two polar groups: the donor group
is a hydrogen atom covalently bound to an electronegative atom (usually O, N, F) while
the acceptor group is another highly electronegative atom.
The criteria used depend on a maximum donor-acceptor distance and a maximum acceptor-donor-H angle.
Hydrogen bonds criteria can be specified.

For more information about chemfiles selection language, please see
http://chemfiles.org/chemfiles/latest/selections.html

Usage:
  cfiles hbonds [options] <trajectory>
  cfiles hbonds (-h | --help)

Examples:
  cfiles hbonds water.xyz --cell 15:15:25 --guess-bonds
  cfiles hbonds in.pdb --donors=="bonds: type(#1) == O and type(#2) == H"
  cfiles hbonds protein.pdb --acceptors=="atoms: type N" --angle 20.0

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
  --donors=<sel>                selection to use for the donors. This must be a
                                selection of size 2 and type bonds The second atom
				should be the hydrogen atom.
                                [default: bonds: type(#2) == H]
  --acceptors=<sel>             selection to use for the acceptors. This must be a
                                selection of size 1.
                                [default: atoms: type O or type N or type F]
  --distance=<dist>             distance criterion to use for the hydrogen bond detection.
                                'dist' is the donor-acceptor maximum distance in angstroms.
                                [default: 3.0]
  --angle=<angle>               angle criterion to use for the hydrogen bond detection.
                                'angle' is the acceptor-donor-hydrogen maximum angle in degrees.
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

    options.selection_acceptor = args.at("--acceptors").asString();
    options.selection_donor = args.at("--donors").asString();

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
        options.distance = string2double(args.at("--distance").asString());
    }	
    if (args.at("--angle")) {
        options.angle = string2double(args.at("--angle").asString())*pi/180;
    }	

    return options;
}


std::string HBonds::description() const {
    return "compute hydrogen bonds network";
}

int HBonds::run(int argc, const char* argv[]) {
    auto options = parse_options(argc, argv);

    selection_donor_ = Selection(options.selection_donor);
    if (selection_donor_.size() != 2) {
        throw CFilesError("Can not use a selection for donors with size different than 2.");
    }

    auto selection_string = split(selection_donor_.string(),':');
    if (selection_string[0] != "bonds") {
        throw CFilesError("'Bonds' type selection compulsory for donors.");
    }

    selection_acceptor_ = Selection(options.selection_acceptor);
    if (selection_acceptor_.size() != 1) {
        throw CFilesError("Can not use a selection for acceptors with size larger than 1.");
    }

    auto infile = Trajectory(options.trajectory, 'r', options.format);
    std::ofstream outfile(options.outfile, std::ios::out);
    if (outfile.is_open()) {
        outfile << "#Hydrogen bond network in trajectory " << options.trajectory << std::endl;
        outfile << "# Selection: acceptors: " << options.selection_acceptor;
        outfile << " and donors: " << options.selection_donor << std::endl;
        outfile << "# Criteria:" << std::endl;
        outfile << "# donor-acceptor distance < " << options.distance << " angstroms" << std::endl;
        outfile << "# acceptor-donor-H angle < " << options.angle * 180 / pi << " degrees" << std::endl;
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
        outfile << "# Acceptor (name index)   Donor (name index)   Hydrogen (name index)  : Distance D-A    Angle A-D-H" << std::endl;

        auto positions = frame.positions();
        auto cell = frame.cell();

        auto matched = selection_donor_.evaluate(frame);
        for (auto match: matched) {
            assert(match.size() == 2);

            size_t donor = static_cast<size_t>(-1);
            size_t hydrogen = static_cast<size_t>(-1);

            if (frame.topology()[match[0]].type() == "H") {
                donor = match[1];
                hydrogen = match[0];
            } else if (frame.topology()[match[1]].type() == "H") {
                donor = match[0];
                hydrogen = match[1];
            } else
            {
                throw CFilesError("Invalid donors selection: there is no hydrogen atom.");
            }

            auto acceptors = selection_acceptor_.list(frame);
            for (auto acceptor: acceptors) {
                if (acceptor != donor && frame.topology()[acceptor].type() != "H") {
                    auto r_da = cell.wrap(positions[acceptor] - positions[donor]);
                    auto distance = norm(r_da);
                    auto r_dh = cell.wrap(positions[hydrogen] - positions[donor]);
                    auto theta = angle(r_da, r_dh);
                    if (distance < options.distance && theta < options.angle) {
                        outfile << frame.topology()[acceptor].name() << " " << acceptor << "   ";
                        outfile << frame.topology()[donor].name() << " " << donor << "   ";
                        outfile << frame.topology()[hydrogen].name() << " " << hydrogen << "  : ";
                        outfile << distance << "    ";
                        outfile << theta * 180 / pi << std::endl;
                    }
                }
            }
        }
    }

    return 0;
}
