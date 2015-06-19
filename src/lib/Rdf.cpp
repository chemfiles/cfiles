/*
 * chrp, an analysis frontend for the Chemharp library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#include <string>

#include "docopt/docopt.h"
#include <Chemharp.hpp>

#include "Rdf.hpp"
#include "Errors.hpp"

std::vector<std::string> split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> elems;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

static const char OPTIONS[] =
R"(chrp rdf: radial distribution function calculations
Usage:
  chrp rdf [options] <trajectory>
  chrp rdf (-h | --help)

Options:
  -h --help                     show this help
  -o <file>, --output=<file>    write result to <file>
  --max=<max>                   maximal distance to use [default: 10]
  --start=<n>                   first step [default: 0]
  --end=<n>                     last step (-1 is the last step) [default: -1]
  --stride=<n>                  use a step every <n> steps [default: 1]
  -b <n>, --bins=<n>            number of bins in the histogram [default: 200]
  -t <path>, --topology=<path>  path to an alternative topology file
  -c <cell>, --cell=<cell>      alternative unit cell. <cell> should be formated
                                using the <a:b:c:alpha:beta:gamma> or <a:b:c> or
                                <L> format. This option set <max> to L/2.
)";

static void parse_options(int argc, char** argv, rdf_options& options) {
    auto args = docopt::docopt(OPTIONS, {argv, argv + argc}, true, "");
    options.infile = args["<trajectory>"].asString();

    if (args["--output"]){
        options.outfile = args["--output"].asString();
    } else {
        options.outfile = options.infile + ".rdf";
    }

    options.rmax = stod(args["--max"].asString());
    options.nbins = stol(args["--bins"].asString());
    options.start = stol(args["--start"].asString());
    options.end = stol(args["--end"].asString());
    options.stride = stol(args["--stride"].asString());

    if (args["--topology"]){
        options.topology = args["--topology"].asString();
    } else {
        options.topology = "";
    }

    if (args["--cell"]) {
		auto cell_strs = split(args["--cell"].asString(), ':');
		if (cell_strs.size()==1) {
			options.cell.push_back(stod(cell_strs[0]));
			options.cell.push_back(stod(cell_strs[0]));
			options.cell.push_back(stod(cell_strs[0]));
		} else if (cell_strs.size()==3) {
			options.cell.push_back(stod(cell_strs[0]));
			options.cell.push_back(stod(cell_strs[1]));
			options.cell.push_back(stod(cell_strs[2]));
		} else if (cell_strs.size()==3) {
			options.cell.push_back(stod(cell_strs[0]));
			options.cell.push_back(stod(cell_strs[1]));
			options.cell.push_back(stod(cell_strs[2]));
            options.cell.push_back(stod(cell_strs[3]));
            options.cell.push_back(stod(cell_strs[4]));
            options.cell.push_back(stod(cell_strs[5]));
		} else {
			throw chrp_exception("The cell option should have 1, 3 or 6 values.");
		}

        double L = std::min(options.cell[0], std::min(options.cell[1], options.cell[2]));
        options.rmax = L/2;
	}
}

std::string Rdf::description() const {
    return "Compute radial distribution functions";
}

std::string Rdf::help() const {
    return OPTIONS;
}

static const double pi = 3.141592653589793238463;

int Rdf::run(int argc, char** argv) {
    parse_options(argc, argv, options_);

    histogram_ = Histogram<double>(options_.nbins, 0, options_.rmax);
    result_ = std::vector<double>(histogram_.size(), 0);

    harp::Trajectory file(options_.infile);

    if (options_.cell.size() == 3) {
        file.cell(harp::UnitCell(options_.cell[0], options_.cell[2], options_.cell[2]));
    } else if (options_.cell.size() == 6) {
        file.cell(harp::UnitCell(
            options_.cell[0], options_.cell[2], options_.cell[2],
            options_.cell[3], options_.cell[4], options_.cell[5]
        ));
    }

    if (options_.topology != "") {
        file.topology(options_.topology);
    }

    size_t start=options_.start, end=options_.end, stride=options_.stride;
    if (end == -1){
        end = file.nsteps();
    }
    if (start > end) {
        chrp_exception("Can not have 'start' bigger than 'end'.");
    }

    for (size_t i=start; i<end; i+=stride) {
        harp::Frame frame = file.read();
        accumulate(frame);
    }
    finish();
    write(options_.outfile);
    return 0;
}

#include <iostream>
using namespace std;

void Rdf::accumulate(harp::Frame& frame) {
    auto positions = frame.positions();
    auto cell = frame.cell();
    size_t npairs = 0;

	for(size_t i=0; i<frame.natoms(); i++){
		// TODO: add selection to only compute RDF on some atoms
		auto ri = positions[i];
        npairs++;
		for(size_t j=i+1; j<frame.natoms(); j++){
			// TODO: add selection to only compute RDF on some atoms

            auto rj = positions[j];
            double rij = norm(cell.wrap(ri - rj));
            if (rij < options_.rmax){
			    histogram_.insert(rij);
            }
		}
	}
	nsteps_++;

    double V = cell.volume();
    double dr = histogram_.bin_size();
    double norm = 1;

    if (V > 0) {
        norm = 2*pi * npairs * npairs / V * dr;
    } else {
        norm = 2*pi * dr;
    }

    histogram_.normalize([&norm, &dr](size_t i, double val){
        double rr = (i + 0.5)*dr * (i + 0.5)*dr;
        return val / (norm * rr);
    });

    for (size_t i=0; i<histogram_.size(); i++){
        result_[i] += histogram_[i];
        histogram_[i] = 0;
    }
}

void Rdf::finish() {
    for (size_t i=0; i<result_.size(); i++){
        result_[i] /= nsteps_;
    }
}

void Rdf::write(const std::string& filename) {
    std::ofstream outfile(filename, std::ios::out);
    if(outfile.is_open()) {
        outfile << "# Radial distribution function for file " << options_.infile << std::endl;

        double dr = histogram_.bin_size();
        for (size_t i=0; i<result_.size(); i++){
            outfile << i * dr << "  " << result_[i] << "\n";
        }

        outfile << std::endl;
    } else {
        throw chrp_exception("Could not open the " + filename + "file.");
    }
}
