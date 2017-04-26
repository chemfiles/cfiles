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
R"(Compute the density profile of particles along a given axis. Selections for the particles
can be specified using the chemfiles selection language. It is possible to provide an alternative
unit cell or topology for the trajectory file if they are not defined in the trajectory format.
The axis can be specified using a coordinate vector (e.g. z axis would be (0,0,1)).

For more information about chemfiles selection language, please see
http://chemfiles.org/chemfiles/latest/selections.html

Usage:
  cfiles density [options] <trajectory>
  cfiles density (-h | --help)

Examples:
  cfiles density water.xyz --cell 15:15:25 --guess-bonds --axis=1,1,1
  cfiles density in.pdb --selection="x>3"

Options:
  -h --help                     show this help
  -o <file>, --output=<file>    write result to <file>. This default to the
                                trajectory file name with the `_dp.dat` extension.
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
                                and <stride> optional. Default is to use all
                                steps from the input; starting at 0, ending at
                                the last step, and with a stride of 1.
  -s <sel>, --selection=<sel>   selection to use for the particles. This must be a
                                selection of size 1.
                                [default: atoms: all]
  --axis=<axis>                 axis along which the density profile will be 
                                computed. It should be either one of 'X','Y','Z'
                                or a vector defining the axis (e.g. 1:1:1). 
                                [default: Z]
  --profile=<str>               type of density profile. Currently implemented are
                                "along_axis" and "radial". [default: along_axis]
  -p <n>, --points=<n>          number of points in the profile [default: 200]
  --max=<m>                     maximum distance in the profile. [default: 10]
  --min=<m>                     minimum distance in the profile.
                                If "--profile=radial", min is 0.
                                If "--profile=along_axis", min is -max if not precised.
)";

Averager<double> DensityProfile::setup(int argc, const char* argv[]) {
    auto options_str = command_header("density", DensityProfile().description()) + "\n";
    options_str += "Laura Scalfi <laura.scalfi@ens.fr>\n";
    options_str += OPTIONS;
    auto args = docopt::docopt(options_str, {argv, argv + argc}, true, "");

    AveCommand::parse_options(args);

    options_.npoints = string2long(args.at("--points").asString());
    options_.selection = args.at("--selection").asString();

    selection_ = Selection(options_.selection);
    if (selection_.size() != 1) {
        throw CFilesError("Can not use a selection with size different than 1.");
    }

    if (args.at("--output")){
        options_.outfile = args.at("--output").asString();
    } else {
        options_.outfile = AveCommand::options().trajectory + "_dp.dat";
    }

    options_.type_profile = args.at("--profile").asString();

    if (args.at("--axis")) {
        auto splitted = split(args.at("--axis").asString(),':');
        if (splitted.size() == 1) {
            axis_ = Axis(args.at("--axis").asString());
        } else if (splitted.size() == 3) {
            auto a = string2double(splitted[0]);
            auto b = string2double(splitted[1]);
            auto c = string2double(splitted[2]);
            options_.axis = vector3d(a,b,c);
            if (options_.axis == vector3d(0,0,0)) {
                throw CFilesError("Null axis does not make sense");
            }
            axis_ = Axis(options_.axis[0], options_.axis[1], options_.axis[2]);
        } else {
            throw CFilesError("Vector should be of size 3");
        }
    } 

    if (args.at("--max")) {
        options_.max = string2double(args.at("--max").asString());
        options_.min = - options_.max;
        if (args.at("--min")) {
            options_.min = string2double(args.at("--min").asString());
        }
        if (options_.type_profile == "radial") {
            options_.min = 0;
        }
        return Averager<double>(options_.npoints, options_.min, options_.max);
    } else {
        throw CFilesError("Please enter a maximum value for the distance values");
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
        if (options_.type_profile == "along_axis") {
            z = axis_.projection(cell.wrap(positions[i]));
        } else if (options_.type_profile == "radial") {
            z = axis_.radial(cell.wrap(positions[i]));
        }
        try {
            profile.insert(z);
        } catch (const OutOfBoundsError& e) {
            std::cout << e.what() << std::endl;
        }
    }
}

void DensityProfile::finish(const Histogram<double>& profile) {
    auto max = *std::max_element(profile.begin(), profile.end());

    if (max == 0) {
        throw CFilesError("No particles found in the '" + selection_.string() + "' selection found.");
    }

    std::ofstream outfile(options_.outfile, std::ios::out);
    if(outfile.is_open()) {
        outfile << "# Density profile in trajectory " << AveCommand::options().trajectory << std::endl;
        outfile << "# along axis " << axis_.coordinates()[0] << ' ' << axis_.coordinates()[1] << ' ';
        outfile << axis_.coordinates()[2] << std::endl;
        outfile << "# Selection: " << options_.selection << std::endl;

        double dr = profile.bin_size();
        for (size_t i=0; i<profile.size(); i++){
            outfile << profile.min() + i * dr << "  " << profile[i] / max << "\n";
        }
    } else {
        throw CFilesError("Could not open the '" + options_.outfile + "' file.");
    }
}
