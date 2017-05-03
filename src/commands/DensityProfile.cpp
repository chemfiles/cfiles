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
language. It is possible to provide an alternative unit cell or topology 
for the trajectory file if they are not defined in the trajectory format.
The axis can be specified using a coordinate vector (e.g. z axis would be 
(0,0,1)).

It is possible to compute 2D profiles. This is possible by specifying 2 axis 
(see --axis and --radial options). Other options (--points, --max, --min, 
--origin) may accept two values, one for each axis. If only one is specified, 
the same value will be used for both axis (see Examples). The output is
a 2D histogram with the first dimension (raws) being the first axis and the 
second dimension (columns) the second axis. If two axis of the same type are
used (e.g. twice --axis option), the order will be the one the user gave. If 
the axis types are different (e.g. --axis and --radial), the --axis will be
first.  

For more information about chemfiles selection language, please see
http://chemfiles.org/chemfiles/latest/selections.html

Usage:
  cfiles density [options] <trajectory>
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
  --origin=<coord>              coordinates for the origin of the axis (used 
                                only for radial profiles).
                                [default: 0:0:0]  
  -p <n>, --points=<n>          number of points in the profile [default: 200]
  --max=<max>                   maximum distance in the profile. [default: 10]
  --min=<min>                   minimum distance in the profile.
                                For radial profiles, <min> must be positive.
                                For linear profiles, <min> is set to -<max> 
                                if not precised.
)";

Axis axis_parse(std::string axis) {
    auto splitted = split(axis,':');
    if (splitted.size() == 1) {
        return Axis(axis);
    } else if (splitted.size() == 3) {
        auto a = string2double(splitted[0]);
        auto b = string2double(splitted[1]);
        auto c = string2double(splitted[2]);
        return Axis(a, b, c);
    } else {
        throw CFilesError("Axis for density profile should be of size 3");
    }
}

Averager<double> DensityProfile::setup(int argc, const char* argv[]) {
    auto options_str = command_header("density", DensityProfile().description()) + "\n";
    options_str += "Laura Scalfi <laura.scalfi@ens.fr>\n";
    options_str += OPTIONS;
    auto args = docopt::docopt(options_str, {argv, argv + argc}, true, "");

    AveCommand::parse_options(args);

    size_t n_axis = 0;

    if (args.at("--axis")) {
        if (args["--axis"].isString()) {
            options_.type_profile[n_axis] = 1;
            n_axis ++;
            axis_x_ = axis_parse(args.at("--axis").asString());
        } else if (args["--axis"].isStringList()) {
            if (args.at("--axis").asStringList().size() == 2) {
                options_.type_profile[n_axis] = 1;
                n_axis ++;
                axis_x_ = axis_parse(args.at("--axis").asStringList()[0]);
                options_.type_profile[n_axis] = 1;
                n_axis ++;
                axis_y_ = axis_parse(args.at("--axis").asStringList()[1]);
            } else {
                throw CFilesError("Too many axis were given");
            }
        }
    }

    if (args.at("--radial")) {
        if (args["--radial"].isString()) {
            options_.type_profile[n_axis] = 2;
            n_axis ++;
            if (n_axis > 2) {
                throw CFilesError("Too many axis were given");
            }
            axis_x_ = axis_parse(args.at("--radial").asString());
        } else if (args["--radial"].isStringList()) {
            if (args.at("--radial").asStringList().size() == 2) {
                options_.type_profile[n_axis] = 2;
                n_axis ++;
                if (n_axis > 2) {
                    throw CFilesError("Too many axis were given");
                }
                axis_x_ = axis_parse(args.at("--axis").asStringList()[0]);
                options_.type_profile[n_axis] = 2;
                n_axis ++;
                if (n_axis > 2) {
                    throw CFilesError("Too many axis were given");
                }
                axis_y_ = axis_parse(args.at("--radial").asStringList()[1]);
            } else {
                throw CFilesError("Too many axis were given");
            }
        }
    }


    options_.npoints = string2long(args.at("--points").asString());
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

    if (args.at("--origin")) {
        auto splitted = split(args.at("--origin").asString(),':');
        if (splitted.size() != 3) {
            throw CFilesError("Origin for density profile should be a vector of size 3");
        }
        auto a = string2double(splitted[0]);
        auto b = string2double(splitted[1]);
        auto c = string2double(splitted[2]);
        options_.origin = vector3d(a,b,c);
    }

    if (args.at("--max")) {
        options_.max = string2double(args.at("--max").asString());
        options_.min = - options_.max;
        if (args.at("--min")) {
            options_.min = string2double(args.at("--min").asString());
        }
        if (options_.type_profile[0] == 2) {
            options_.min = 0;
        }
        return Averager<double>(options_.npoints, options_.min, options_.max);
    } else {
        throw CFilesError("Please enter a maximum value for the distances in density profile");
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
        double z;
        if (options_.type_profile[0] == 1) {
            z = axis_x_.projection(cell.wrap(positions[i]));
        } else if (options_.type_profile[0] == 2) {
            z = axis_x_.radial(cell.wrap(positions[i]-options_.origin));
        }
        try {
            profile.insert(z);
        } catch (const OutOfBoundsError& e) {
            std::cout << e.what() << std::endl;
        }
    }
}

void DensityProfile::finish(const Histogram<double>& profile) {
    std::ofstream outfile(options_.outfile, std::ios::out);
    auto axis = axis_x_.get_coordinates();
    if(outfile.is_open()) {
        outfile << "# Density profile in trajectory " << AveCommand::options().trajectory << std::endl;
        outfile << "# along axis " << axis[0] << ' ' << axis[1] << ' ' << axis[2] << std::endl;
        outfile << "# Selection: " << options_.selection << std::endl;

        double dr = profile.dx();
        for (size_t i=0; i<profile.size(); i++){
            if (options_.type_profile[0] == 1) {
                outfile << profile.min_x() + i * dr << "  " << profile[i] << "\n";
            } else if (options_.type_profile[0] == 2) {
                auto r = profile.min_x() + (i + 0.5) * dr;
                if (r == 0) {
                    r = dr / 1000; // use a small r compared to dr to avoid Nan in the output
                }
                outfile << profile.min_x() + i * dr << "  " << profile[i] / r << "\n";
            }
        }
    } else {
        throw CFilesError("Could not open the '" + options_.outfile + "' file.");
    }
}
