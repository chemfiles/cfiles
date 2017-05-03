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

    n_axis_ = 0;

    if (args.at("--axis")) {
        if (args.at("--axis").asStringList().size() == 0) {
        } else if (args.at("--axis").asStringList().size() == 1) {
            n_axis_ ++;
            axis_x_ = axis_parse(args.at("--axis").asStringList()[0]);
            options_.type_profile[n_axis_ - 1] = 1;
        } else if (args.at("--axis").asStringList().size() == 2) {
            n_axis_ ++;
            options_.type_profile[n_axis_ - 1] = 1;
            axis_x_ = axis_parse(args.at("--axis").asStringList()[0]);
            n_axis_ ++;
            options_.type_profile[n_axis_ - 1] = 1;
            axis_y_ = axis_parse(args.at("--axis").asStringList()[1]);
        } else {
            throw CFilesError("Too many axis were given");
        }
    }

    if (args.at("--radial")) {
        if (args.at("--radial").asStringList().size() == 0) {
        } else if (args.at("--radial").asStringList().size() == 1) {
            n_axis_ ++;
            if (n_axis_ == 1) {
                axis_x_ = axis_parse(args.at("--radial").asStringList()[0]);
            } else if (n_axis_ == 2) {
                axis_y_ = axis_parse(args.at("--radial").asStringList()[0]);
            } else { 
                throw CFilesError("Too many axis were given");
            }
            options_.type_profile[n_axis_ - 1] = 2;
        } else if (args.at("--radial").asStringList().size() == 2) {
            if (n_axis_ == 2) {
                throw CFilesError("Too many axis were given");
            } else {
                n_axis_ ++;
                options_.type_profile[n_axis_ - 1] = 2;
                axis_x_ = axis_parse(args.at("--radial").asStringList()[0]);
                n_axis_ ++;
                options_.type_profile[n_axis_ - 1] = 2;
                axis_y_ = axis_parse(args.at("--radial").asStringList()[0]);
            }
        }
    }

    if (args.at("--points")) {
        auto splitted = split(args.at("--points").asString(),':');
        if (splitted.size() == 1) {
            options_.npoints[0] = string2long(splitted[0]);
            options_.npoints[1] = string2long(splitted[0]);
        } else if (splitted.size() == 2) {
            if (n_axis_ < 2) {
                throw CFilesError("More --points options than axis");
            }
            options_.npoints[0] = string2long(splitted[0]);
            options_.npoints[1] = string2long(splitted[1]);
        } else {
            throw CFilesError("Too many arguments for --points option");
        }
    }

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
        if (splitted.size() != 3 and splitted.size() != 6) {
            throw CFilesError("Origin for density profile should be one or two vectors of size 3");
        }
        if (splitted.size() == 3) {
            auto a = string2double(splitted[0]);
            auto b = string2double(splitted[1]);
            auto c = string2double(splitted[2]);
            options_.origin[0] = vector3d(a, b, c);
            options_.origin[1] = vector3d(a, b, c);
        } else if (splitted.size() == 6) {
            if (options_.type_profile[0] != 2 or options_.type_profile[1] != 2) {
                throw CFilesError("More --origin options than radial axis");
            }
            auto a1 = string2double(splitted[0]);
            auto b1 = string2double(splitted[1]);
            auto c1 = string2double(splitted[2]);
            auto a2 = string2double(splitted[3]);
            auto b2 = string2double(splitted[4]);
            auto c2 = string2double(splitted[5]);
            options_.origin[0] = vector3d(a1, b1, c1);
            options_.origin[1] = vector3d(a2, b2, c2);
        }
    }

    if (args.at("--max")) {
        auto splitted = split(args.at("--max").asString(),':');
        if (splitted.size() == 1) {
            options_.max[0] = string2double(splitted[0]);
            options_.max[1] = string2double(splitted[0]);
        } else if (splitted.size() == 2) {
            if (n_axis_ < 2) {
                throw CFilesError("More --max options than axis");
            }
            options_.max[0] = string2double(splitted[0]);
            options_.max[1] = string2double(splitted[1]);
        } else {
            throw CFilesError("Too many arguments for --points option");
        }
        options_.min[0] = - options_.max[0];
        options_.min[1] = - options_.max[1];

        if (args.at("--min")) {
             auto splitted = split(args.at("--min").asString(),':');
             if (splitted.size() == 1) {
                 options_.min[0] = string2double(splitted[0]);
                 options_.min[1] = string2double(splitted[0]);
             } else if (splitted.size() == 2) {
                 if (n_axis_ < 2) {
                     throw CFilesError("More --min options than axis");
                 }
                 options_.min[0] = string2double(splitted[0]);
                 options_.min[1] = string2double(splitted[1]);
             } else {
                 throw CFilesError("Too many arguments for --points option");
             }
        }
        if (options_.type_profile[0] == 2) {
            options_.min[0] = 0;
        } else if (options_.type_profile[1] == 2) {
            options_.min[1] = 0;
        }
    } else {
        throw CFilesError("Please enter a maximum value for the distances in density profile");
    }
        
    if (n_axis_ == 1) {
        return Averager<double>(options_.npoints[0], options_.min[0], options_.max[0]);
    } else if (n_axis_ == 2) {
        return Averager<double>(options_.npoints[0], options_.min[0], options_.max[0], options_.npoints[1], options_.min[1], options_.max[1]);
    } else {
        throw CFilesError("No axis or too many axis have been defined");
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
        double x;
        double y;
        if (options_.type_profile[0] == 1) {
            x = axis_x_.projection(cell.wrap(positions[i]));
        } else if (options_.type_profile[0] == 2) {
            x = axis_x_.radial(cell.wrap(positions[i]-options_.origin[0]));
        }
        if (options_.type_profile[1] == 1) {
            y = axis_y_.projection(cell.wrap(positions[i]));
        } else if (options_.type_profile[1] == 2) {
            y = axis_y_.radial(cell.wrap(positions[i]-options_.origin[1]));
        }
        try {
            if (n_axis_ == 1) {
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
    auto axis_x = axis_x_.get_coordinates();
    auto axis_y = axis_y_.get_coordinates();
    if (outfile.is_open()) {
        outfile << "# Density profile in trajectory " << AveCommand::options().trajectory << std::endl;
        outfile << "# along axis " << axis_x[0] << ' ' << axis_x[1] << ' ' << axis_x[2] << std::endl;
        if (n_axis_ == 2) {
            outfile << "# and along axis " << axis_y[0] << ' ' << axis_y[1] << ' ' << axis_y[2] << std::endl;
        }
        outfile << "# Selection: " << options_.selection << std::endl;

        double dx = profile.dx();
        double dy = profile.dy();
        if (n_axis_ == 1) {
            for (size_t i=0; i<profile.size(); i++){
                if (options_.type_profile[0] == 1) {
                    outfile << profile.min_x() + i * dx << "  " << profile[i] << "\n";
                } else if (options_.type_profile[0] == 2) {
                    auto r = profile.min_x() + (i + 0.5) * dx;
                    if (r == 0) {
                        r = dx / 1000; // use a small r compared to dr to avoid Nan in the output
                    }
                    outfile << profile.min_x() + i * dx << "  " << profile[i] / r << "\n";
                }
            }
        } else {
            outfile << "# X Y Density" << std::endl;
            size_t nx = profile.n_x();
            size_t ny = profile.n_y();
            size_t min_x = profile.min_x();
            size_t min_y = profile.min_y();

            for (size_t i=0; i<nx; i++){
                for (size_t j=0; j<ny; j++){
                    outfile << min_x + i * dx << "\t" << min_y + j * dy << "\t";
                    if (options_.type_profile[0] == 1 and options_.type_profile[1] == 1) {
                        outfile << profile(i,j) << "\n";
                    } else if (options_.type_profile[0] == 2 xor options_.type_profile[1] == 2) {
                            auto r = min_y + (j + 0.5) * dy;
                            if (r == 0) {
                                r = dy / 1000; // use a small r compared to dr to avoid Nan in the output
                            }
                        outfile << profile(i,j) / r << "\n";
                    } else if (options_.type_profile[0] == 2 and options_.type_profile[1] == 2) {
                            auto rx = min_x + (i + 0.5) * dx;
                            if (rx == 0) {
                                rx = dx / 1000; // use a small r compared to dr to avoid Nan in the output
                            }
                            auto ry = min_y + (j + 0.5) * dy;
                            if (ry == 0) {
                                ry = dy / 1000; // use a small r compared to dr to avoid Nan in the output
                            }
                        outfile << profile(i,j) / (rx * ry);
                    }
                }
            }
        }
    } else {
        throw CFilesError("Could not open the '" + options_.outfile + "' file.");
    }
}
