// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <docopt/docopt.h>
#include <unordered_set>
#include <fstream>

#include "Rdf.hpp"
#include "Errors.hpp"
#include "utils.hpp"
#include "warnings.hpp"

using namespace chemfiles;
constexpr double PI = 3.141592653589793238463;

/// Get the radius of the biggest inscribed sphere in the unit cell
static double biggest_sphere_radius(const UnitCell& cell);

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
  --center=<sel/positions>      compute rdf with respect to a single center
                                point instead of using pair distances. The
                                center can either be a fixed position given as
                                a 3D vector (0:2:0), or a selection. If a
                                selection is used, the center of mass of
                                selected atoms will be used as center point.
  --max=<max>                   maximal distance to use. If a custom unit cell
                                is present (--cell option) and this option is
                                not, the radius of the biggest inscribed sphere
                                is used as maximal distance [default: 10]
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

    if (args["--center"]) {
        options_.center = args["--center"].asString();
    }

    options_.rmax = string2double(args["--max"].asString());
    options_.npoints = string2long(args["--points"].asString());
    options_.selection = args["--selection"].asString();

    auto begin = argv;
    auto end = argv + argc;
    auto max_options = std::find_if(begin, end, [](const std::string& arg){
        return arg.substr(0, 5) == "--max";
    });
    if (AveCommand::options().custom_cell && max_options == end) {
        options_.rmax = biggest_sphere_radius(AveCommand::options().cell);
    }

    selection_ = Selection(options_.selection);
    if (selection_.size() > 2) {
        throw CFilesError("Can not use a selection with more than two atoms in RDF.");
    }

    if (!options_.center.empty()) {
        auto center = split(options_.center, ':');
        if (center.size() == 3) {
            center_ = Vector3D(
                string2double(center[0]),
                string2double(center[1]),
                string2double(center[2])
            );
        } else {
            center_sel_ = Selection(options_.center);
            if (selection_.size() != 1) {
                throw CFilesError("Can not use a selection with more than one atoms with a center.");
            }
        }
    }

    coordination_ = Averager<double>(options_.npoints, 0, options_.rmax);
    return Averager<double>(options_.npoints, 0, options_.rmax);
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
    check_rmax(frame);

    auto positions = frame.positions();
    auto cell = frame.cell();
    size_t n_first = 0;
    size_t n_second = 0;

    bool use_center = false;
    auto center = Vector3D();
    if (center_) {
        // The center point is given as a vector
        use_center = true;
        center = *center_;
    } else if (center_sel_) {
        // The center point is the center of mass of a selection
        use_center = true;
        auto& positions = frame.positions();

        auto mass = 0.0;
        for (auto i: center_sel_->list(frame)) {
            mass += frame[i].mass();
            // TODO: Should this use PBC?
            center = center + frame[i].mass() * positions[i];
        }
        center = center / mass;
    }


    if (selection_.size() == 1) {
        // Use the same selection for both atoms in the pair
        auto matched = selection_.list(frame);
        n_first = matched.size();

        if (use_center) {
            // The center point is given as a vector
            auto& cell = frame.cell();
            auto& positions = frame.positions();
            n_second = 1;
            for (auto i: matched) {
                auto rij = center - positions[i];
                cell.wrap(rij);
                auto d = rij.norm();
                if (d < options_.rmax){
                    histogram.insert(d);
                }
            }
        } else {
            // Use the same selection for both atoms in the pair
            n_second = matched.size();
            for (auto i: matched) {
                for (auto j: matched) {
                    if (i == j) continue;

                    auto rij = frame.distance(i, j);
                    if (rij < options_.rmax){
                        histogram.insert(rij);
                    }
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

            auto rij = frame.distance(i, j);
            if (rij < options_.rmax){
                histogram.insert(rij);
            }
        }

        n_first = first_particles.size();
        n_second = second_particles.size();
    }

    if (n_first == 0 || n_second == 0) {
        warn_once(
            "No pair corresponding to '" + selection_.string() + "' found."
        );
        return;
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

    if (use_center) {
        factor = 4 * PI * n_first / volume;
    } else {
        double rho = (n_first + n_second) / volume;
        double alpha = static_cast<double>(n_second) / static_cast<double>(n_first + n_second);
        factor = alpha * 4 * PI * rho;
    }
    for (size_t i=1; i<histogram.size(); i++){
        auto r = (i + 0.5) * dr;
        coordination_[i] = coordination_[i - 1] + factor * histogram[i] * r * r * dr;
    }
    coordination_.step();
}

void Rdf::check_rmax(const chemfiles::Frame& frame) const {
    auto r_sphere = biggest_sphere_radius(frame.cell());
    if (r_sphere < options_.rmax) {
        warn_once(fmt::format(
            "The maximal distance (--max option) is too big for this cell.\n"
            "The cell contains values up to {:.2f} and the max distance is {}.",
            r_sphere, options_.rmax
        ));
    }
}

double biggest_sphere_radius(const UnitCell& cell) {
    auto matrix = cell.matrix();
    auto a = Vector3D(matrix[0][0], matrix[1][0], matrix[2][0]);
    auto b = Vector3D(matrix[0][1], matrix[1][1], matrix[2][1]);
    auto c = Vector3D(matrix[0][2], matrix[1][2], matrix[2][2]);
    // Make sure we have an upper triangular matrix
    assert(matrix[1][0] == 0);
    assert(matrix[2][0] == 0);
    assert(matrix[2][1] == 0);

    // Plan normal vectors
    auto na = cross(b, c);
    auto nb = cross(c, a);
    auto nc = cross(a, b);

    auto ra = std::abs(dot(na, a)) / (2 * na.norm());
    auto rb = std::abs(dot(nb, b)) / (2 * nb.norm());
    auto rc = std::abs(dot(nc, c)) / (2 * nc.norm());

    return std::min(ra, std::min(rb, rc));
}
