/* cfiles, an analysis frontend for the Chemfiles library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#include "Convert.hpp"
#include "utils.hpp"

#include "docopt/docopt.h"
#include "chemfiles.hpp"
using namespace chemfiles;

static const char OPTIONS[] =
R"(cfiles convert: trajectory conversion

Convert trajectories from one format to another, automatically guessing the format to use
based on the files extension. The '--input-format' and '--output-format' can be used to
force the format.

Usage:
  cfiles convert [options] <input> <output>
  cfiles convert (-h | --help)

Options:
  -h --help                     show this help
  --input-format=<format>       force the input file format to be <format>
  --output-format=<format>      force the output file format to be <format>
  -t <path>, --topology=<path>  path to an alternative topology file for the input
  -c <cell>, --cell=<cell>      use alternative unit cell for the input. <cell>
                                should be formated using the <a:b:c:alpha:beta:gamma>
                                or <a:b:c> or <L> format.
)";

struct convert_options {
    std::string in_file;
    std::string out_file;
    std::string in_format;
    std::string out_format;
    std::vector<double> cell;
    std::string topology;
};

static convert_options parse_options(int argc, const char* argv[]) {
    convert_options options;
    auto args = docopt::docopt(OPTIONS, {argv, argv + argc}, true, "");

    options.in_file = args["<input>"].asString();
    options.out_file = args["<output>"].asString();

    if (args["--input-format"]){
        options.in_format = args["--input-format"].asString();
    } else {
        options.in_format = "";
    }

    if (args["--output-format"]){
        options.out_format = args["--output-format"].asString();
    } else {
        options.out_format = "";
    }

    if (args["--topology"]){
        options.topology = args["--topology"].asString();
    } else {
        options.topology = "";
    }

    if (args["--cell"]) {
		options.cell = parse_cell(args["--cell"].asString());
	}
    return options;
}


std::string Convert::description() const {
    return "Convert trajectories between formats";
}

std::string Convert::help() const {
    return OPTIONS;
}

int Convert::run(int argc, const char* argv[]) {
    auto options = parse_options(argc, argv);

    auto infile = Trajectory(options.in_file, 'r', options.in_format);
    auto outfile = Trajectory(options.out_file, 'w', options.out_format);

    if (options.cell.size() == 3) {
        infile.set_cell(UnitCell(options.cell[0], options.cell[2], options.cell[2]));
    } else if (options.cell.size() == 6) {
        infile.set_cell(UnitCell(
            options.cell[0], options.cell[2], options.cell[2],
            options.cell[3], options.cell[4], options.cell[5]
        ));
    }

    if (options.topology != "") {
        infile.set_topology(options.topology);
    }

    for (size_t i=0; i<infile.nsteps(); i++) {
        auto frame = infile.read();
        outfile.write(frame);
    }

    return 0;
}
