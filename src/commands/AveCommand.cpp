// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

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
  -c <cell>, --cell=<cell>      alternative unit cell. <cell> should be formated
                                using one of the <a:b:c:α:β:γ> or <a:b:c> or <L>
                                formats. This option set <max> to L/2.
  --start=<n>                   first step [default: 0]
  --end=<n>                     last step (-1 for the input end) [default: -1]
  --stride=<n>                  use a step every <n> steps [default: 1])";

void AveCommand::parse_options(const std::map<std::string, docopt::value>& args) {
    options_.trajectory = args.at("<trajectory>").asString();
    options_.start = stol(args.at("--start").asString());
    options_.end = stol(args.at("--end").asString());
    options_.stride = stol(args.at("--stride").asString());
    options_.guess_bonds = args.at("--guess-bonds").asBool();

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

    size_t start = options_.start;
    size_t end = options_.end;
    size_t stride = options_.stride;
    if (end == static_cast<size_t>(-1)) {
        end = file.nsteps();
    }
    if (start > end) {
        throw CFilesError("Can not have 'start' bigger than 'end'.");
    }

    for (size_t i=start; i<end; i+=stride) {
        auto frame = file.read();
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
