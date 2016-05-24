/* cfiles, an analysis frontend for the Chemfiles library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

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
                                single selection ("name O") or two selections
                                separated by a comma ("name O, name H")
                                [default: all]
  -t <path>, --topology=<path>  path to an alternative topology file
  -c <cell>, --cell=<cell>      alternative unit cell. <cell> should be formated
                                using the <a:b:c:alpha:beta:gamma> or <a:b:c> or
                                <L> format. This option set <max> to L/2.
  --start=<n>                   first step [default: 0]
  --end=<n>                     last step (-1 is the last step) [default: -1]
  --stride=<n>                  use a step every <n> steps [default: 1]
  --max=<max>                   maximal distance to use [default: 10]
  -p <n>, --points=<n>          number of points in the histogram [default: 200]

)";

static rdf_options parse_options(int argc, const char* argv[]) {
    auto args = docopt::docopt(OPTIONS, {argv, argv + argc}, true, "");

    rdf_options options;
    options.infile = args["<trajectory>"].asString();

    if (args["--output"]){
        options.outfile = args["--output"].asString();
    } else {
        options.outfile = options.infile + ".rdf";
    }

    options.rmax = stod(args["--max"].asString());
    options.npoints = stol(args["--points"].asString());
    options.start = stol(args["--start"].asString());
    options.end = stol(args["--end"].asString());
    options.stride = stol(args["--stride"].asString());

    if (args["--selection"]){
        auto selections = split(args["--selection"].asString(), ',');
        if (selections.size() == 1) {
            options.selection_i = selections[0];
            options.selection_j = selections[0];
        } else if (selections.size() == 2) {
            options.selection_i = selections[0];
            options.selection_j = selections[1];
        } else {
            throw CFilesError("Wrong size for selection string. Only two selection are allowed.");
        }
    } else {
        options.selection_i = "all";
        options.selection_j = "all";
    }

    if (args["--topology"]){
        options.topology = args["--topology"].asString();
    } else {
        options.topology = "";
    }

    if (args["--cell"]) {
		options.cell = parse_cell(args["--cell"].asString());
        double L = std::min(options.cell[0], std::min(options.cell[1], options.cell[2]));
        options.rmax = L/2;
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

    auto file = Trajectory(options_.infile);

    if (options_.cell.size() == 3) {
        file.set_cell(UnitCell(options_.cell[0], options_.cell[2], options_.cell[2]));
    } else if (options_.cell.size() == 6) {
        file.set_cell(UnitCell(
            options_.cell[0], options_.cell[2], options_.cell[2],
            options_.cell[3], options_.cell[4], options_.cell[5]
        ));
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

    sel_i = Selection(options_.selection_i);
    sel_j = Selection(options_.selection_j);

    for (size_t i=start; i<end; i+=stride) {
        auto frame = file.read();
        accumulate(frame);
    }
    finish();
    write(options_.outfile);
    return 0;
}

void Rdf::accumulate(Frame& frame) {
    auto positions = frame.positions();
    auto cell = frame.cell();
    size_t npairs = 0;

    std::vector<Bool> matched_i, matched_j;
    if (options_.selection_i == options_.selection_j) {
        // Only evaluate the selection once
        matched_i = sel_i.evaluate(frame);
        matched_j = matched_i;
    } else {
        matched_i = sel_i.evaluate(frame);
        matched_j = sel_j.evaluate(frame);
    }

	for(size_t i=0; i<frame.natoms(); i++){
	    if (!matched_i[i]) continue;
		auto& ri = positions[i];

		for(size_t j=0; j<frame.natoms(); j++){
            if (i == j) continue;
			if (!matched_j[j]) continue;

            auto& rj = positions[j];
            double rij = norm(cell.wrap(ri - rj));
            if (rij < options_.rmax){
                npairs++;
			    histogram_.insert_at(rij);
            }
		}
	}

    double V = cell.volume();
    if (V > 0) V = 1;

    double dr = histogram_.bin_size();
    double norm = 8 * pi * npairs / V * dr;

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
        outfile << "# Radial distribution function for file " << options_.infile << std::endl;
        outfile << "# Atoms: " << options_.selection_i << " & " << options_.selection_j << std::endl;

        double dr = histogram_.bin_size();
        for (size_t i=0; i<result_.size(); i++){
            outfile << i * dr << "  " << result_[i] << "\n";
        }
        outfile << std::endl;
    } else {
        throw CFilesError("Could not open the '" + filename + "' file.");
    }
}
