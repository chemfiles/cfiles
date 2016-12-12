// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#include <docopt/docopt.h>
#include <algorithm>

#include "Angles.hpp"
#include "Errors.hpp"
#include "geometry.hpp"

using namespace chemfiles;

static const double pi = 3.141592653589793238463;
static const std::string OPTIONS =
R"(cfiles angles: angles histograms

Compute angles distribution along a trajectory. The angle to study can be
specified using the chemfiles selection language. It is also possible to provide
an alternative topology or unit cell when this information is not present in the
trajectory.

Usage:
  cfiles angles [options] <trajectory>
  cfiles angles (-h | --help)

Options:
  -h --help                     show this help
  -o <file>, --output=<file>    write result to <file>. This default to the
                                trajectory file name with the `.ang` extension.
  -s <sel>, --selection=<sel>   selection to use for the atoms. This must be a
                                selection of size 3 (for angles) or 4 (for
                                dihedral angles) [default: angles: all]
  -p <n>, --points=<n>          number of points in the histogram [default: 200])";


std::string AngleDistribution::description() const {
    return "Compute angles distribution";
}

std::string AngleDistribution::help() const {
    return OPTIONS;
}

void AngleDistribution::setup(int argc, const char* argv[], Histogram<double>& histogram) {
    auto options = std::string(OPTIONS) + AveCommand::AVERAGE_OPTIONS;
    auto args = docopt::docopt(options, {argv, argv + argc}, true, "");

    AveCommand::parse_options(args);

    if (args["--output"]){
        options_.outfile = args["--output"].asString();
    } else {
        options_.outfile = AveCommand::options().trajectory + ".ang";
    }

    options_.npoints = stol(args["--points"].asString());
    options_.selection = args["--selection"].asString();

    selection_ = Selection(options_.selection);
    if (selection_.size() == 3) {
        histogram = Histogram<double>(options_.npoints, 0, pi);
    } else if (selection_.size() == 4) {
        histogram = Histogram<double>(options_.npoints, 0, 2*pi);
    } else {
        throw CFilesError("Can not use a selection with less than three atoms in angle distribution.");
    }
}

void AngleDistribution::finish(const Histogram<double>& histogram) {
    auto max = *std::max_element(histogram.begin(), histogram.end());

    std::ofstream outfile(options_.outfile, std::ios::out);
    if(outfile.is_open()) {
        outfile << "# Angles distribution in trajectory " << AveCommand::options().trajectory << std::endl;
        outfile << "# Selection: " << options_.selection << std::endl;

        double dr = histogram.bin_size();
        for (size_t i=0; i<histogram.size(); i++){
            outfile << i * dr << "  " << histogram[i] / max << "\n";
        }
    } else {
        throw CFilesError("Could not open the '" + options_.outfile + "' file.");
    }
}

void AngleDistribution::accumulate(const Frame& frame, Histogram<double>& histogram) {
    auto positions = frame.positions();
    auto cell = frame.cell();

    auto matched = selection_.evaluate(frame);
    for (auto match: matched) {
        assert(match.size() == 3 || match.size() == 4);

        if (match.size() == 3) {
            auto i = match[0];
            auto j = match[1];
            auto k = match[2];
            auto r21 = cell.wrap(positions[i] - positions[j]);
            auto r23 = cell.wrap(positions[k] - positions[j]);
            auto theta = angle(r21, r23);
            histogram.insert(theta);
        } else if (match.size() == 4) {
            auto i = match[0];
            auto j = match[1];
            auto k = match[2];
            auto m = match[3];
            auto r12 = cell.wrap(positions[j] - positions[i]);
            auto r23 = cell.wrap(positions[k] - positions[j]);
            auto r34 = cell.wrap(positions[m] - positions[k]);
            auto phi = dihedral(r12, r23, r34);
            histogram.insert(phi);
        }
    }
}
