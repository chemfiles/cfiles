// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <docopt/docopt.h>
#include <unordered_set>
#include <fstream>

#include "Rdf.hpp"
#include "Errors.hpp"
#include "utils.hpp"

using namespace chemfiles;

static const double PI = 3.141592653589793238463;
static const char OPTIONS[] =
R"(Compute radial pair distribution function (often denoted g(r)) and running
coordination number. The pair of particles to use can be specified using the
chemfiles selection language. It is possible to provide an alternative unit
cell or topology for the trajectory file if they are not defined in the
trajectory format.

For more information about chemfiles selection language, please see
http://chemfiles.github.io/chemfiles/latest/selections.html

Usage:
  cfiles rdf [options] <trajectory>
  cfiles rdf (-h | --help)

Examples:
  cfiles rdf water.tng -s "name O" --max=8.5 --output=rdf-O-O.dat
  cfiles rdf butane.tng -s "pairs: name(#1) C and name(#2) H"
  cfiles rdf methane.xyz --cell 15:15:25 --guess-bonds --points=150
  cfiles rdf result.xtc --topology=initial.mol --topology-format=PDB
  cfiles rdf simulation.pdb --steps=10000::100 -o partial-rdf.dat

Options:
  -h --help                     show this help
  -o <file>, --output=<file>    write result to <file>. This default to the
                                trajectory file name with the `.rdf.dat`
                                extension.
  -s <sel>, --selection=<sel>   selection to use for the atoms. This can be a
                                single selection ("name O") or a selection of
                                two atoms ("pairs: name(#1) O and name(#2) H")
                                [default: all]
  --max=<max>                   maximal distance to use [default: 10]
  -p <n>, --points=<n>          number of points in the histogram [default: 200])";

std::string Rdf::description() const {
    return "compute radial distribution functions";
}

Averager<double> Rdf::setup(int argc, const char* argv[]) {
    auto options = command_header("rdf", Rdf().description());
    options += "Guillaume Fraux <guillaume@fraux.fr>\n\n";
    options += std::string(OPTIONS) + AveCommand::AVERAGE_OPTIONS;
    auto args = docopt::docopt(options, {argv, argv + argc}, true, "");
    AveCommand::parse_options(args);

    if (args["--output"]){
        options_.outfile = args["--output"].asString();
    } else {
        options_.outfile = AveCommand::options().trajectory + ".rdf.dat";
    }

    options_.rmax = string2double(args["--max"].asString());
    options_.npoints = string2long(args["--points"].asString());
    options_.selection = args["--selection"].asString();

    if (AveCommand::options().custom_cell) {
        auto& cell = AveCommand::options().cell;
        double L = std::min(cell.a(), std::min(cell.b(), cell.c()));
        options_.rmax = L/2;
	}

    selection_ = Selection(options_.selection);
    if (selection_.size() > 2) {
        throw CFilesError("Can not use a selection with more than two atoms in RDF.");
    } else {
        coordination_ = Averager<double>(options_.npoints, 0, options_.rmax);
        return Averager<double>(options_.npoints, 0, options_.rmax);
    }
}

void Rdf::finish(const Histogram<double>& histogram) {
    coordination_.average();

    std::ofstream outfile(options_.outfile, std::ios::out);
    if(!outfile.is_open()) {
        throw CFilesError("Could not open the '" + options_.outfile + "' file.");
    }

    outfile << "# Radial distribution function in trajectory " << AveCommand::options().trajectory << std::endl;
    outfile << "# Using selection: " << options_.selection << std::endl;
    outfile << "# r\tg(r)\tN(r) " << std::endl;

    for (size_t i=0; i<histogram.size(); i++){
        outfile << histogram.first_coord(i) << "\t" << histogram[i] << "\t" << coordination_[i] << "\n";
    }
}

void Rdf::accumulate(const Frame& frame, Histogram<double>& histogram) {
    auto positions = frame.positions();
    auto cell = frame.cell();
    size_t n_first = 0;
    size_t n_second = 0;

    if (selection_.size() == 1) {
        // If we have a single atom selection, use it for both atoms of the
        // pairs
        auto matched = selection_.list(frame);
        n_first = matched.size();
        n_second = matched.size();
        for (auto i: matched) {
            for (auto j: matched) {
                if (i == j) continue;

                auto rij = norm(cell.wrap(positions[j] - positions[i]));
                if (rij < options_.rmax){
                    histogram.insert(rij);
                }
            }
    	}
    } else {
        // If we have a pair selection, use it directly
        assert(selection_.size() == 2);
        auto matched = selection_.evaluate(frame);
        std::unordered_set<size_t> first_particles;
        std::unordered_set<size_t> second_particles;

        for (auto match: matched) {
    		auto i = match[0];
            auto j = match[1];

            first_particles.insert(i);
            second_particles.insert(j);

            auto rij = norm(cell.wrap(positions[j] - positions[i]));
            if (rij < options_.rmax){
                histogram.insert(rij);
            }
    	}

        n_first = first_particles.size();
        n_second = second_particles.size();
    }

    if (n_first == 0 || n_second == 0) {
        throw CFilesError(
            "No pair corresponding to the '" + selection_.string() + "' selection found."
        );
    }

    // Normalize the rdf to be 1 at long distances
    double volume = cell.volume();
    if (volume <= 0) {volume = 1;}

    double dr = histogram.first().dr;
    double factor = n_first * n_second / volume;

    histogram.normalize([factor, dr](size_t i, double val){
        double r = (i + 0.5) * dr;
        return val / (4 * PI * factor * dr * r * r);
    });

    double rho = (n_first + n_second) / volume;
    for (size_t i=1; i<histogram.size(); i++){
        auto r = (i + 0.5) * dr;
        coordination_[i] = coordination_[i - 1] + 4 * PI * rho * histogram[i] * r * r * dr;
    }
    coordination_.step();
}
