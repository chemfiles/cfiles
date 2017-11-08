// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <iostream>
#include <sstream>

#include <docopt/docopt.h>
#include <chemfiles.hpp>

#include "Info.hpp"
#include "Errors.hpp"
#include "utils.hpp"

using namespace chemfiles;

static const char OPTIONS[] =
R"(Get various information and metadata from a trajectory.

Usage:
  cfiles info [options] <input>
  cfiles info (-h | --help)

Examples:
    cfiles info water.xyz
    cfiles info --guess-bonds --step 4 water.xyz

Options:
  -h --help                     show this help
  --guess-bonds                 guess the bonds in the input
  --step=<step>                 give informations about the frame at <step>
                                [default: 0]
)";

static Info::Options parse_options(int argc, const char* argv[]) {
    auto options_str = command_header("info", Info().description());
    options_str += "Guillaume Fraux <guillaume@fraux.fr>\n\n";
    options_str += OPTIONS;
    auto args = docopt::docopt(options_str, {argv, argv + argc}, true, "");

    Info::Options options;
    options.input = args["<input>"].asString();
    options.guess_bonds = args.at("--guess-bonds").asBool();
    auto step = args.at("--step").asLong();

    if (step >= 0) {
        options.step = static_cast<size_t>(step);
    } else {
        throw CFilesError("step must be positive");
    }

    return options;
}

std::string Info::description() const {
    return "get information on a trajectory";
}

int Info::run(int argc, const char* argv[]) {
    auto options = parse_options(argc, argv);
    auto input = Trajectory(options.input, 'r', "");

    std::stringstream output;

    output << "information for " << options.input << std::endl;

    output << std::endl << "global:" << std::endl;
    output << "    steps = " << input.nsteps() <<  std::endl;

    if (input.nsteps() > options.step) {
        auto frame = input.read_step(options.step);
        output << std::endl << "frame " << frame.step() << ":" << std::endl;
        output << "    atoms = " << frame.size() <<  std::endl;

        if (options.guess_bonds) {
            frame.guess_topology();
        }

        auto& topology = frame.topology();
        output << "    bonds = " << topology.bonds().size() <<  std::endl;
        output << "    angles = " << topology.angles().size() <<  std::endl;
        output << "    dihedrals = " << topology.dihedrals().size() <<  std::endl;
        output << "    residues = " << topology.residues().size() <<  std::endl;
    }


    std::cout << output.str();

    return 0;
}
