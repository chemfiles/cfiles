// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <docopt/docopt.h>
#include <sstream>

#include "AveCommand.hpp"
#include "Errors.hpp"
#include "utils.hpp"

using namespace chemfiles;

const std::string AveCommand::AVERAGE_OPTIONS = R"(
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
                                the number of steps for <end> and 1 for <stride>.)";

void AveCommand::parse_options(const std::map<std::string, docopt::value>& args) {
    options_.trajectory = args.at("<trajectory>").asString();
    options_.guess_bonds = args.at("--guess-bonds").asBool();

    if (args.at("--steps")) {
        options_.steps = steps_range::parse(args.at("--steps").asString());
    }

    if (args.at("--topology")) {
        if (options_.guess_bonds) {
            throw CFilesError("Can not use both '--topology' and '--guess-bonds'");
        }
        options_.topology = args.at("--topology").asString();
    }

    if (args.at("--format")) {
        options_.format = args.at("--format").asString();
    }

    if (args.at("--topology-format")) {
        if (options_.topology == "") {
            throw CFilesError("Can not use '--topology-format' without a '--topology'");
        }
        options_.topology_format = args.at("--topology-format").asString();
    }

    if (args.at("--cell")) {
        options_.custom_cell = true;
		options_.cell = parse_cell(args.at("--cell").asString());
	}
}

int AveCommand::run(int argc, const char* argv[]) {
    histogram_ = setup(argc, argv);

    auto file = Trajectory(options_.trajectory);
    if (options_.custom_cell) {
        file.set_cell(options_.cell);
    }
    if (options_.topology != "") {
        file.set_topology(options_.topology, options_.topology_format);
    }

    for (auto step: options_.steps) {
        if (step >= file.nsteps()) {
            break;
        }
        auto frame = file.read_step(step);
        if (options_.guess_bonds) {
            frame.guess_topology();
        }
        accumulate(frame, histogram_);
        histogram_.step();
    }
    histogram_.average();

    finish(histogram_);
    return 0;
}
