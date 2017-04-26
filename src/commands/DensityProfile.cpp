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
  --selection=<sel>             selection to use for the particles. This must be a
                                selection of size 1.
                                [default: atoms: all]
  --axis=<axis>                 axis along which the density profile will be computed.
                                To choose from 'x', 'y' or 'z'. [default: z]
  -p <n>, --points=<n>          number of points in the profile [default: 200]
  --max=<m>                     maximum distance in the profile. The default is the box
                                length in the chosen direction.
                                Be careful if your box size changes
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

    if (args.at("--axis")) {
        options_.axis_str = args.at("--axis").asString();
    }
    axis_ = Axis(options_.axis_str);

    if (args.at("--max")) {
        return Averager<double>(options_.npoints, 0, options_.max);
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
        
        double z = axis_.projection(cell.wrap(positions[i]));

        try {
            profile.insert(z);
        } catch (const std::string e) {
            std::cout << e << std::endl;
        }
    }
}

void DensityProfile::finish(const Histogram<double>& profile) {
    auto max = *std::max_element(profile.begin(), profile.end());

    if (max == 0) {
        throw CFilesError("No particles found in the '" + selection_.string() + "' selection found.");
    }

    char ax[3] = {'x','y','z'};

    std::ofstream outfile(options_.outfile, std::ios::out);
    if(outfile.is_open()) {
        outfile << "# Density profile in trajectory " << AveCommand::options().trajectory << std::endl;
        outfile << "# along axis " << options_.axis_str << ',';
        outfile << "# Selection: " << options_.selection << std::endl;

        double dr = profile.bin_size();
        for (size_t i=0; i<profile.size(); i++){
            outfile << i * dr << "  " << profile[i] / max << "\n";
        }
    } else {
        throw CFilesError("Could not open the '" + options_.outfile + "' file.");
    }
}
