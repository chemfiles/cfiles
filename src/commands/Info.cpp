// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <iostream>
#include <sstream>

#include <fmt/format.h>
#include <fmt/ostream.h>
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
  --format=<format>             force the input file format to be <format>
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

    if (args.at("--format")) {
        options.format = args.at("--format").asString();
    }

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
    auto input = Trajectory(options.input, 'r', options.format);

    std::stringstream output;
    fmt::print(output, "file = {}\n", options.input);
    fmt::print(output, "steps = {}\n", input.nsteps());

    if (input.nsteps() > options.step) {
        auto frame = input.read_step(options.step);
        fmt::print(output, "\n[frame(step={})]\n", frame.step());

        auto& cell = frame.cell();
        auto lengths = cell.lengths();
        auto angles = cell.angles();
        fmt::print(output, "cell = [({}, {}, {}), ({}, {}, {})]\n",
            lengths[0], lengths[1], lengths[2],
            angles[0], angles[1], angles[2]
        );

        if (options.guess_bonds) {
            frame.guess_bonds();
        }

        auto& topology = frame.topology();
        fmt::print(output, "atoms_count = {}\n", frame.size());
        fmt::print(output, "bonds_count = {}\n", topology.bonds().size());
        fmt::print(output, "angles_count = {}\n", topology.angles().size());
        fmt::print(output, "dihedrals_count = {}\n", topology.dihedrals().size());
        fmt::print(output, "impropers_count = {}\n", topology.impropers().size());
        fmt::print(output, "residues_count = {}\n", topology.residues().size());
    }

    std::cout << output.str();

    return 0;
}
