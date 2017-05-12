// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <docopt/docopt.h>
#include <sstream>
#include <fstream>
#include <iostream>

#include "DensityProfile.hpp"
#include "Errors.hpp"
#include "utils.hpp"
#include "geometry.hpp"

using namespace chemfiles;

static const double pi = 3.141592653589793238463;
static const char OPTIONS[] =
R"(Compute the density profile of particles along a given axis or radially. 
The output for the radial density profile is normalized by r.
Selections for the particles can be specified using the chemfiles selection 
language. It is possible to provide an alternative unit cell or topology for the
trajectory file if they are not defined in the trajectory format. The axis can 
be specified using a coordinate vector (e.g. z axis would be (0,0,1)).

It is also possible to compute 2D profiles by specifying 2 axis (see --axis and
--radial options). Other options (--points, --max, --min) may accept two values,
one for each axis. If only one is specified, the same value will be used for 
both axis (see Examples). The output is a 2D histogram with the first dimension
being the first axis and the second dimension the second axis. If two axis of 
the same type are used (e.g. twice --axis option), the order will be the one the
user gave. If the axis types are different (e.g. --axis and --radial), the 
--axis will be first. Two axis of type radial are forbidden.  

For more information about chemfiles selection language, please see
http://chemfiles.org/chemfiles/latest/selections.html

Usage:
  cfiles density [options] <trajectory> [--axis=<axis>...] [--radial=<axis>...]
  cfiles density (-h | --help)

Examples:
  cfiles density water.xyz --cell 15:15:25 --guess-bonds --axis=1:1:1
  cfiles density in.pdb --selection="x > 3" --points=500
  cfiles density nt.pdb --radial=Z --max=3 --origin=0:0:2
  cfiles density nt.pdb --profile=Z --radial=Z --max=10:5 --origin=0:0:2

Options:
  -h --help                     show this help
  -o <file>, --output=<file>    write result to <file>. This default to the
                                trajectory file name with the `.density.dat` 
                                extension.
  --format=<format>             force the input file format to be <format>
  -t <path>, --topology=<path>  alternative topology file for the input
  --topology-format=<format>    use <format> as format for the topology file
  --guess-bonds                 guess the bonds in the input
  -c <cell>, --cell=<cell>      alternative unit cell. <cell> format is one of
                                <a:b:c:α:β:γ> or <a:b:c> or <a>. 'a', 'b' and
                                'c' are in angstroms, 'α', 'β', and 'γ' are in
                                degrees.
  --steps=<steps>               steps to use from the input. <steps> format
                                is <start>:<end>[:<stride>] with <start>, <end>
                                and <stride> optional. The used steps goes from
                                <start> to <end> (excluded) by steps of
                                <stride>. The default values are 0 for <start>,
                                the number of steps for <end> and 1 for 
                                <stride>.
  -s <sel>, --selection=<sel>   selection to use for the particles. This must
                                be a selection of size 1. [default: atoms: all]
  --axis=<axis>...              computes a linear density profile along <axis>.
                                It should be either one of 'X','Y','Z'
                                or a vector defining the axis (e.g. 1:1:1). 
  --radial=<axis>...            computes a radial density profile using the 
                                distance to <axis>.
                                It should be either one of 'X','Y','Z'
                                or a vector defining the axis (e.g. 1:1:1). 
  --origin=<coord>              coordinates for the origin of the axis (only 
                                relevant for radial profiles). [default: 0:0:0]
  -p <n>, --points=<n>          number of points in the profile [default: 200]
  --max=<max>                   maximum distance in the profile. [default: 10]
  --min=<min>                   minimum distance in the profile. [default: 0]
                                For radial profiles, <min> must be positive.
)";

Averager<double> DensityProfile::setup(int argc, const char* argv[]) {
    auto options_str = command_header("density", DensityProfile().description()) + "\n";
    options_str += "Laura Scalfi <laura.scalfi@ens.fr>\n";
    options_str += OPTIONS;
    auto args = docopt::docopt(options_str, {argv, argv + argc}, true, "");

    AveCommand::parse_options(args);

    options_.selection = args.at("--selection").asString();
    selection_ = Selection(options_.selection);
    if (selection_.size() != 1) {
        throw CFilesError("Can not use a selection with size different than 1.");
    }

    if (args.at("--output")){
        options_.outfile = args.at("--output").asString();
    } else {
        options_.outfile = AveCommand::options().trajectory + ".density.dat";
    }

    if (args.at("--axis")) {
        for (auto axis: args.at("--axis").asStringList()) {
            axis_.push_back(Axis::parse(axis, Axis::Linear));
        }
    }

    if (args.at("--radial")) {
        for (auto axis: args.at("--radial").asStringList()) {
            axis_.push_back(Axis::parse(axis, Axis::Radial));
        }
    }

    size_t dimension = dimensionality();

    if (dimension == 0 or dimension > 2) {
        throw CFilesError("No axis or too many axis were given");
    }

    if (args.at("--points")) {
        auto splitted = split(args.at("--points").asString(),':');
        if (splitted.size() == 1) {
            options_.npoints[0] = string2long(splitted[0]);
            options_.npoints[1] = string2long(splitted[0]);
        } else if (splitted.size() == 2) {
            if (dimension < 2) {
                throw CFilesError("More --points options than axis");
            }
            options_.npoints[0] = string2long(splitted[0]);
            options_.npoints[1] = string2long(splitted[1]);
        } else {
            throw CFilesError("Too many arguments for --points option");
        }
    }

    if (args.at("--origin")) {
        auto splitted = split(args.at("--origin").asString(),':');
        if (splitted.size() != 3) {
            throw CFilesError("Origin for density profile should be a vector of size 3");
        }
        auto a = string2double(splitted[0]);
        auto b = string2double(splitted[1]);
        auto c = string2double(splitted[2]);
        options_.origin = vector3d(a, b, c);
    }

    if (args.at("--max")) {
        auto splitted = split(args.at("--max").asString(),':');
        if (splitted.size() == 1) {
            options_.max[0] = string2double(splitted[0]);
            options_.max[1] = string2double(splitted[0]);
        } else if (splitted.size() == 2) {
            if (dimension < 2) {
                throw CFilesError("More --max options than axis");
            }
            options_.max[0] = string2double(splitted[0]);
            options_.max[1] = string2double(splitted[1]);
        } else {
            throw CFilesError("Too many arguments for --max option");
        }
    }

    if (args.at("--min")) {
         auto splitted = split(args.at("--min").asString(),':');
         if (splitted.size() == 1) {
             options_.min[0] = string2double(splitted[0]);
             options_.min[1] = string2double(splitted[0]);
         } else if (splitted.size() == 2) {
             if (dimension < 2) {
                 throw CFilesError("More --min options than axis");
             }
             options_.min[0] = string2double(splitted[0]);
             options_.min[1] = string2double(splitted[1]);
         } else {
             throw CFilesError("Too many arguments for --min option");
         }
    }

    if (options_.min[0] > options_.max[0]) {
        throw CFilesError("Min > Max for first dimension");
    }

    if (options_.min[1] > options_.max[1]) {
        throw CFilesError("Min > Max for second dimension");
    }

    if (dimension == 1) {
        if (axis_[0].is_radial()) {
            if (options_.min[0] < 0) {
                throw CFilesError("Min value for radial axis should be positive");
            }
        }
        return Averager<double>(options_.npoints[0], options_.min[0], options_.max[0]);
    } else {
        assert(dimension == 2);
        if (axis_[0].is_radial()) {
            if (options_.min[1] < 0) {
                throw CFilesError("Min value for radial axis should be positive");
            }
        }
        return Averager<double>(options_.npoints[0], options_.min[0], options_.max[0], options_.npoints[1], options_.min[1], options_.max[1]);
    }
}


std::string DensityProfile::description() const {
    return "compute density profiles";
}

void DensityProfile::accumulate(const chemfiles::Frame& frame, Histogram<double>& profile) {
    auto positions = frame.positions();
    auto cell = frame.cell();

    auto matched = selection_.evaluate(frame);
    for (auto match: matched) {
        assert(match.size() == 1);

        auto i = match[0];
        double x = 0;
        double y = 0;
        if (axis_[0].is_linear()) {
            x = axis_[0].projection(cell.wrap(positions[i]));
        } else {
            assert(axis_[0].is_radial());
            x = axis_[0].radial(cell.wrap(positions[i] - options_.origin));
        }
        if (dimensionality() == 2) {
            if (axis_[1].is_linear()) {
                y = axis_[1].projection(cell.wrap(positions[i]));
            } else {
                assert(axis_[1].is_radial());
                y = axis_[1].radial(cell.wrap(positions[i] - options_.origin));
            }
        }
        try {
            if (dimensionality() == 1) {
                profile.insert(x);
            } else {
                profile.insert(x,y);
            }
        } catch (const OutOfBoundsError& e) {
            std::cout << e.what() << std::endl;
        }
    }
}

void DensityProfile::finish(const Histogram<double>& profile) {
    std::ofstream outfile(options_.outfile, std::ios::out);
    auto axis1 = axis_[0].get_coordinates();
    auto axis2 = axis_[1].get_coordinates();
    if (outfile.is_open()) {
        outfile << "# Density profile in trajectory " << AveCommand::options().trajectory << std::endl;
        outfile << "# along axis " << axis1[0] << ' ' << axis1[1] << ' ' << axis1[2] << std::endl;
        if (dimensionality() == 2) {
            outfile << "# and along axis " << axis2[0] << ' ' << axis2[1] << ' ' << axis2[2] << std::endl;
        }
        outfile << "# Selection: " << options_.selection << std::endl;

        if (dimensionality() == 1) {
            for (size_t i = 0; i < profile.size(); i++){
                if (axis_[0].is_linear()) {
                    outfile << profile.first_coord(i) << "  " << profile[i] << "\n";
                } else {
                    assert(axis_[0].is_radial());
                    outfile << profile.first_coord(i) << "  " << profile[i] / profile.first_coord(i) << "\n";
                }
            }
        } else {
            outfile << "# FirstDimension SecondDimension Density" << std::endl;
            
            for (size_t i = 0; i < profile.first().nbins; i++){
                for (size_t j = 0; j < profile.second().nbins; j++){
                    outfile << profile.first_coord(i) << "\t" << profile.second_coord(j) << "\t";
                    if (axis_[0].is_linear() and axis_[1].is_linear()) {
                        outfile << profile(i,j) << "\n";
                    } else {
                        assert(axis_[0].is_linear() and axis_[1].is_radial());
                        outfile << profile(i,j) / profile.second_coord(j) << "\n";
                    }
                }
            }
        }
    } else {
        throw CFilesError("Could not open the '" + options_.outfile + "' file.");
    }
}
