// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <docopt/docopt.h>
#include <algorithm>
#include <fstream>

#include "Angles.hpp"
#include "Errors.hpp"
#include "warnings.hpp"

using namespace chemfiles;
constexpr double PI = 3.141592653589793238463;

static inline double rad2deg(double value) {
    return value * 180 / PI;
}

static const std::string OPTIONS =
R"(Compute distribution of angles or dihedral angles along a trajectory. The
angle can be specified using the chemfiles selection language. It is possible
to provide an alternative unit cell or topology for the trajectory file if they
are not defined in the trajectory format.

For more information about chemfiles selection language, please see
http://chemfiles.github.io/chemfiles/latest/selections.html

Usage:
  cfiles angles [options] <trajectory>
  cfiles angles (-h | --help)

Examples:
  cfiles angles water.tng -s "angles: name(#1) H and name(#2) O and name(#3) H"
  cfiles angles butane.tng -s "dihedrals: name(#2) C and name(#3) C"
  cfiles angles methane.xyz --cell 15:15:25 --guess-bonds --points=150
  cfiles angles result.xtc --topology=initial.mol --topology-format=PDB
  cfiles angles simulation.pdb --steps=:1000:5 -o partial-angles.dat

Options:
  -h --help                     show this help
  -o <file>, --output=<file>    write result to <file>. This default to the
                                trajectory file name with the `.angles.dat`
                                extension.
  -s <sel>, --selection=<sel>   selection to use for the atoms. This must be a
                                selection of size 3 (for angles) or 4 (for
                                dihedral angles) [default: angles: all]
  -p <n>, --points=<n>          number of points in the histogram [default: 200])";


std::string Angles::description() const {
    return "compute angles and dihedral angles distribution";
}

Averager Angles::setup(int argc, const char* argv[]) {
    auto options = command_header("angles", Angles().description());
    options += "Guillaume Fraux <guillaume@fraux.fr>\n\n";
    options += std::string(OPTIONS) + AveCommand::AVERAGE_OPTIONS;
    auto args = docopt::docopt(options, {argv, argv + argc}, true, "");

    AveCommand::parse_options(args);

    if (args["--output"]) {
        options_.outfile = args["--output"].asString();
    } else {
        options_.outfile = AveCommand::options().trajectory + ".angles.dat";
    }

    options_.npoints = string2long(args["--points"].asString());
    options_.selection = args["--selection"].asString();

    selection_ = Selection(options_.selection);
    if (selection_.size() == 3) {
        return Averager(options_.npoints, 0, PI);
    } else if (selection_.size() == 4) {
        return Averager(options_.npoints, -PI, PI);
    } else {
        throw CFilesError("Can not use a selection with less than three atoms in angle distribution.");
    }
}

void Angles::finish(const Histogram& histogram) {
    double sum = 0;
    for (size_t i=0; i<histogram.size(); i++) {
        sum += rad2deg(histogram.first().width) * histogram[i];
    }

    std::ofstream outfile(options_.outfile, std::ios::out);
    if(outfile.is_open()) {
        outfile << "# Angles distribution in trajectory " << AveCommand::options().trajectory << std::endl;
        outfile << "# Selection: " << options_.selection << std::endl;

        for (size_t i=0; i<histogram.size(); i++) {
            outfile << rad2deg(histogram.first().coord(i)) << "  " << histogram[i] / sum << "\n";
        }
    } else {
        throw CFilesError("Could not open the '" + options_.outfile + "' file.");
    }
}

void Angles::accumulate(const Frame& frame, Histogram& histogram) {
    auto matched = selection_.evaluate(frame);
    if (matched.empty()) {
        warn_once(
            "No angle corresponding to '" + selection_.string() + "' found."
        );
    }

    for (auto match: matched) {
        assert(match.size() == 3 || match.size() == 4);

        if (match.size() == 3) {
            auto theta = frame.angle(match[0], match[1], match[2]);
            histogram.insert(theta);
        } else if (match.size() == 4) {
            auto phi = frame.dihedral(match[0], match[1], match[2], match[3]);
            histogram.insert(phi);
        }
    }
}
