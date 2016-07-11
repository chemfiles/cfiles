/* cfiles, an analysis frontend for the Chemfiles library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
#include <docopt/docopt.h>

#include "AverageCommand.hpp"
#include "Errors.hpp"
#include "utils.hpp"

using namespace chemfiles;

const std::string AverageCommand::AVERAGE_OPTIONS = R"(
  -t <path>, --topology=<path>  alternative topology file for the input
  -c <cell>, --cell=<cell>      alternative unit cell. <cell> should be formated
                                using one of the <a:b:c:α:β:γ> or <a:b:c> or <L>
                                formats. This option set <max> to L/2.
  --start=<n>                   first step [default: 0]
  --end=<n>                     last step (-1 for the input end) [default: -1]
  --stride=<n>                  use a step every <n> steps [default: 1])";

void AverageCommand::parse_options(const std::map<std::string, docopt::value>& args) {
    options_.trajectory = args.at("<trajectory>").asString();
    options_.start = stol(args.at("--start").asString());
    options_.end = stol(args.at("--end").asString());
    options_.stride = stol(args.at("--stride").asString());

    if (args.at("--topology")){
        options_.topology = args.at("--topology").asString();
    } else {
        options_.topology = "";
    }

    if (args.at("--cell")) {
        options_.custom_cell = true;
		options_.cell = parse_cell(args.at("--cell").asString());
	} else {
        options_.custom_cell = false;
    }
}

int AverageCommand::run(int argc, const char* argv[]) {
    setup(argc, argv, histogram_);
    result_ = std::vector<double>(histogram_.size(), 0);

    auto file = Trajectory(options_.trajectory);

    if (options_.custom_cell) {
        file.set_cell(options_.cell);
    }

    if (options_.topology != "") {
        file.set_topology(options_.topology);
    }

    size_t start=options_.start, end=options_.end, stride=options_.stride;
    if (end == static_cast<size_t>(-1)) {
        end = file.nsteps();
    }
    if (start > end) {
        throw CFilesError("Can not have 'start' bigger than 'end'.");
    }

    for (size_t i=start; i<end; i+=stride) {
        nsteps_ += 1;

        // Accumulate intermediate results in results
        auto frame = file.read();
        accumulate(frame, histogram_);

        for (size_t i=0; i<histogram_.size(); i++){
            result_[i] += histogram_[i];
            histogram_[i] = 0;
        }
    }

    for (size_t i=0; i<result_.size(); i++){
        histogram_[i] = result_[i] / nsteps_;
    }

    finish(histogram_);
    return 0;
}