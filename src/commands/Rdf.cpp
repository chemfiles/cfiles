/* cfiles, an analysis frontend for the Chemfiles library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#include <iostream>

#include <string>

#include "docopt/docopt.h"
#include "chemfiles.hpp"
using namespace chemfiles;

#include "Rdf.hpp"
#include "Errors.hpp"
#include "utils.hpp"

static const char OPTIONS[] =
R"(cfiles rdf: compute radial distribution function

Compute pair radial distrubution function (often called g(r)). The pairs of
particles to use can be specified using the chemfiles selection language. It
is possible to provide an alternative topology or unit cell when this
information is not present in the trajectory.

Usage:
  cfiles rdf [options] <trajectory>
  cfiles rdf (-h | --help)

Options:
  -h --help                     show this help
  -o <file>, --output=<file>    write result to <file>. This default to the
                                trajectory file name with the `.rdf` extension.
  -s <sel>, --selection=<sel>   selection to use for the atoms. This can be a
                                single selection ("name O") or a selection of
                                two atoms ("pairs: name($1) O and name($2) H")
                                [default: all]
  -t <path>, --topology=<path>  alternative topology file for the input
  -c <cell>, --cell=<cell>      alternative unit cell. <cell> should be formated
                                using one of the <a:b:c:α:β:γ> or <a:b:c> or <L>
                                formats. This option set <max> to L/2.
  --start=<n>                   first step [default: 0]
  --end=<n>                     last step (-1 for the input end) [default: -1]
  --stride=<n>                  use a step every <n> steps [default: 1]
  --max=<max>                   maximal distance to use [default: 10]
  -p <n>, --points=<n>          number of points in the histogram [default: 200]

)";

static rdf_options parse_options(int argc, const char* argv[]) {
    auto args = docopt::docopt(OPTIONS, {argv, argv + argc}, true, "");

    rdf_options options;
    options.trajectory = args["<trajectory>"].asString();

    if (args["--output"]){
        options.outfile = args["--output"].asString();
    } else {
        options.outfile = options.trajectory + ".rdf";
    }

    options.rmax = stod(args["--max"].asString());
    options.npoints = stol(args["--points"].asString());
    options.start = stol(args["--start"].asString());
    options.end = stol(args["--end"].asString());
    options.stride = stol(args["--stride"].asString());
    options.selection = args["--selection"].asString();

    if (args["--topology"]){
        options.topology = args["--topology"].asString();
    } else {
        options.topology = "";
    }

    if (args["--cell"]) {
        options.custom_cell = true;
		options.cell = parse_cell(args["--cell"].asString());
        double L = std::min(options.cell.a(), std::min(options.cell.b(), options.cell.c()));
        options.rmax = L/2;
	} else {
        options.custom_cell = false;
    }

    return options;
}

std::string Rdf::description() const {
    return "Compute radial distribution functions";
}

std::string Rdf::help() const {
    return OPTIONS;
}

static const double pi = 3.141592653589793238463;

int Rdf::run(int argc, const char* argv[]) {
    options_ = parse_options(argc, argv);

    histogram_ = Histogram<double>(options_.npoints, 0, options_.rmax);
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

    selection_ = Selection(options_.selection);
    if (selection_.size() > 2) {
        throw CFilesError("Can not use a selection with more than two atoms in RDF.");
    }

    for (size_t i=start; i<end; i+=stride) {
        auto frame = file.read();
        accumulate(frame);
    }
    finish();
    write(options_.outfile);
    return 0;
}

void Rdf::accumulate(const Frame& frame) {
    auto positions = frame.positions();
    auto cell = frame.cell();
    size_t npairs = 0;

    if (selection_.size() == 1) {
        // If we have a single atom selection, use it for both atoms of the
        // pairs
        auto matched = selection_.list(frame);
        for (auto i: matched) {
            for (auto j: matched) {
                if (i == j) continue;

                auto rij = norm(cell.wrap(positions[j] - positions[i]));
                if (rij < options_.rmax){
                    histogram_.insert_at(rij);
                    npairs++;
                }
            }
    	}
    } else {
        // If we have a pair selection, use it directly
        assert(selection_.size() == 2);
        auto matched = selection_.evaluate(frame);

        for (auto match: matched) {
    		auto i = match[0];
    		auto j = match[1];
            auto rij = norm(cell.wrap(positions[j] - positions[i]));
            if (rij < options_.rmax){
                histogram_.insert_at(rij);
                npairs++;
            }
    	}
    }

    // Normalize the rdf to be 1 at long distances
    double V = cell.volume();
    if (V > 0) V = 1;

    double dr = histogram_.bin_size();
    double rho = frame.natoms() / V;
    double norm = 1e-6 * 2 * 4 * pi * rho * npairs * dr;

    histogram_.normalize([norm, dr](size_t i, double val){
        double r = (i + 0.5) * dr;
        return val / (norm * r * r);
    });

    for (size_t i=0; i<histogram_.size(); i++){
        result_[i] += histogram_[i];
        histogram_[i] = 0;
    }
    nsteps_++;
}

void Rdf::finish() {
    for (size_t i=0; i<result_.size(); i++){
        result_[i] /= nsteps_;
    }
}

void Rdf::write(const std::string& filename) {
    std::ofstream outfile(filename, std::ios::out);
    if(outfile.is_open()) {
        outfile << "# Radial distribution function in trajectory " << options_.trajectory << std::endl;
        outfile << "# Selection: " << options_.selection << std::endl;

        double dr = histogram_.bin_size();
        for (size_t i=0; i<result_.size(); i++){
            outfile << i * dr << "  " << result_[i] << "\n";
        }
    } else {
        throw CFilesError("Could not open the '" + filename + "' file.");
    }
}
