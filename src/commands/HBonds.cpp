// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <docopt/docopt.h>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include "HBonds.hpp"
#include "Autocorrelation.hpp"
#include "Histogram.hpp"
#include "Errors.hpp"
#include "utils.hpp"
#include "warnings.hpp"

using namespace chemfiles;
constexpr double PI = 3.141592653589793238463;

static const char OPTIONS[] =
R"(Compute list of hydrogen bonds along a trajectory. Selections for the acceptor
and donor atoms can be specified using the chemfiles selection language. It is
possible to provide an alternative unit cell or topology for the trajectory file
if they are not defined in the trajectory format. Hydrogen bonds are defined as
electrostatic attraction between two polar groups: the donor group is a hydrogen
atom covalently bound to an electronegative atom (usually O, N, F) while the
acceptor group is another highly electronegative atom. The criteria used depend
on a maximum donor-acceptor distance and a maximum acceptor-donor-H angle.
Hydrogen bonds criteria can be specified.

For more information about chemfiles selection language, please see
http://chemfiles.org/chemfiles/latest/selections.html

Usage:
  cfiles hbonds [options] <trajectory>
  cfiles hbonds (-h | --help)

Examples:
  cfiles hbonds water.xyz --cell 15:15:25 --guess-bonds
  cfiles hbonds in.pdb --donors="bonds: type(#1) == O and type(#2) == H"
  cfiles hbonds protein.pdb --acceptors="atoms: type N" --angle 20.0

Options:
  -h --help                     show this help
  -o <file>, --output=<file>    write result to <file>. This default to the
                                trajectory file name with the `.hbonds.dat`
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
  --donors=<sel>                selection to use for the donors. This must be a
                                selection of size 2, with the hydrogen atom as
                                second atom. [default: bonds: type(#2) == H]
  --acceptors=<sel>             selection to use for the acceptors. This must
                                be a selection of size 1.
                                [default: atoms: type O or type N or type F]
  --distance=<distance>         distance criterion to use for the hydrogen bond
                                detection. <distance> is the donor-acceptor
                                maximum distance in angstroms. [default: 3.5]
  --angle=<angle>               angle criterion to use for the hydrogen bond
                                detection. <angle> is the acceptor-donor-hydrogen
                                maximum angle in degrees. [default: 30.0]
  --histogram=<output>          accumulate the hydrogen bond histogram as a
                                function of (r, theta) and output it to the
                                given <ouput> file.
  -p <n>, --points=<n>          number of points in the histogram [default: 200]
  --autocorrelation=<output>    compute the hydrogen bond existence
                                autocorrelation and output it to the given
                                <ouput> file. This can be used to retrieve the
                                lifetime of hydrogen bonds.
)";

struct hbond {
    size_t donor;
    size_t hydrogen;
    size_t acceptor;
};

bool operator==(const hbond& lhs, const hbond& rhs) {
    return (lhs.donor == rhs.donor && lhs.hydrogen == rhs.hydrogen && lhs.acceptor == rhs.acceptor);
}

inline void hash_combine(size_t& hash, size_t value) {
    hash ^= std::hash<size_t>()(value) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
}

namespace std {
    template <> struct hash<hbond> {
        size_t operator()(const hbond& bond) const {
            size_t hash = 0;
            hash_combine(hash, bond.donor);
            hash_combine(hash, bond.hydrogen);
            hash_combine(hash, bond.acceptor);
            return hash;
        }
    };
}

static HBonds::Options parse_options(int argc, const char* argv[]) {
    auto options_str = command_header("hbonds", HBonds().description()) + "\n";
    options_str += "Laura Scalfi <laura.scalfi@ens.fr>\n";
    options_str += OPTIONS;
    auto args = docopt::docopt(options_str, {argv, argv + argc}, true, "");

    HBonds::Options options;
    options.trajectory = args.at("<trajectory>").asString();
    options.guess_bonds = args.at("--guess-bonds").asBool();

    options.acceptor_selection = args.at("--acceptors").asString();
    options.donor_selection = args.at("--donors").asString();

    options.distance = string2double(args.at("--distance").asString());
    options.angle = string2double(args.at("--angle").asString()) * PI / 180;
    options.npoints = string2long(args["--points"].asString());

    if (args.at("--output")) {
        options.outfile = args.at("--output").asString();
    } else {
        options.outfile = options.trajectory + ".hbonds.dat";
    }

    if (args.at("--autocorrelation")) {
        options.autocorr_output = args.at("--autocorrelation").asString();
        options.autocorrelation = true;
    } else {
        options.autocorrelation = false;
    }

    if (args.at("--histogram")) {
        options.histogram_output = args.at("--histogram").asString();
        options.histogram = true;
    } else {
        options.histogram = false;
    }

    if (args.at("--steps")) {
        options.steps = steps_range::parse(args.at("--steps").asString());
    }

    if (args.at("--format")) {
        options.format = args.at("--format").asString();
    }

    if (args.at("--topology")) {
        if (options.guess_bonds) {
            throw CFilesError("Can not use both '--topology' and '--guess-bonds'");
        }
        options.topology = args.at("--topology").asString();
    }

    if (args.at("--topology-format")) {
        if (options.topology == "") {
            throw CFilesError("Can not use '--topology-format' without a '--topology'");
        }
        options.topology_format = args.at("--topology-format").asString();
    }

    if (args.at("--cell")) {
        options.custom_cell = true;
        options.cell = parse_cell(args.at("--cell").asString());
    }

    return options;
}

std::string HBonds::description() const {
    return "compute hydrogen bonds using distance/angle criteria";
}

int HBonds::run(int argc, const char* argv[]) {
    auto options = parse_options(argc, argv);

    auto donors = Selection(options.donor_selection);
    if (donors.size() != 2) {
        throw CFilesError("Can not use a selection for donors with size that is not 2.");
    }

    auto acceptors = Selection(options.acceptor_selection);
    if (acceptors.size() != 1) {
        throw CFilesError("Can not use a selection for acceptors with size larger than 1.");
    }

    std::ofstream outfile(options.outfile, std::ios::out);
    if (!outfile.is_open()) {
        throw CFilesError("Could not open the '" + options.outfile + "' file.");
    }
    fmt::print(outfile, "# Hydrogen bonds in {}\n", options.trajectory);
    fmt::print(outfile, "# Between '{}' and '{}'\n", options.acceptor_selection, options.donor_selection);

    auto infile = Trajectory(options.trajectory, 'r', options.format);
    if (options.custom_cell) {
        infile.set_cell(options.cell);
    }

    if (options.topology != "") {
        infile.set_topology(options.topology, options.topology_format);
    }

    auto histogram = Histogram(options.npoints, 0, options.distance, options.npoints, 0, options.angle * 180 / PI);
    auto existing_bonds = std::unordered_map<hbond, std::vector<float>>();
    size_t used_steps = 0;
    for (auto step: options.steps) {
        if (step >= infile.nsteps()) {
            break;
        }
        auto frame = infile.read_step(step);
        if (options.guess_bonds) {
            frame.guess_bonds();
        }

        auto bonds = std::unordered_set<hbond>();
        auto matched = donors.evaluate(frame);
        if (matched.empty()) {
            warn("no atom matching the donnor selection at step " + std::to_string(step));
        }

        for (auto match: matched) {
            assert(match.size() == 2);

            size_t donor = match[0];
            size_t hydrogen = match[1];

            if (frame[hydrogen].type() != "H") {
                warn_once(
                    "the second atom in the donors selection might not be an "
                    "hydrogen (expected type H, got type " + frame[hydrogen].type() + ")"
                );
            }

            auto acceptors_list = acceptors.list(frame);
            if (acceptors_list.empty()) {
                warn("no atom matching the acceptor selection at step " + std::to_string(step));
            }

            for (auto acceptor: acceptors_list) {
                if (acceptor != donor && frame.topology()[acceptor].type() != "H") {
                    auto distance = frame.distance(acceptor, donor);
                    auto theta = frame.angle(acceptor, donor, hydrogen);
                    if (distance < options.distance && theta < options.angle) {
                        bonds.emplace(hbond{donor, hydrogen, acceptor});
                        if (options.histogram) {
                            histogram.insert(distance, theta * 180 / PI);
                        }
                    }
                }
            }
        }

        fmt::print(outfile, "# step n_bonds\n");
        fmt::print(outfile, "{} {}\n", step, bonds.size());
        fmt::print(outfile, "# Donnor Hydrogen Acceptor\n", step);
        for (auto& bond: bonds) {
            fmt::print(outfile, "{} {} {}\n", bond.donor, bond.hydrogen, bond.acceptor);
        }

        if (options.histogram) {
            auto max = *std::max_element(histogram.begin(), histogram.end());
            histogram.normalize([max](size_t, double value) {
                return value / max;
            });

            std::ofstream outhist(options.histogram_output, std::ios::out);
            if (!outhist.is_open()) {
                throw CFilesError("Could not open the '" + options.histogram_output + "' file.");
            }

            fmt::print(outhist, "# Hydrogen bonds density histogram in {}\n", options.trajectory);
            fmt::print(outhist, "# Between '{}' and '{}'\n", options.acceptor_selection, options.donor_selection);
            fmt::print(outhist, "# r theta density\n", options.acceptor_selection, options.donor_selection);

            for (size_t i = 0; i < histogram.first().nbins; i++){
                for (size_t j = 0; j < histogram.second().nbins; j++){
                    fmt::print(
                        outhist,
                        "{} {} {}\n",
                        histogram.first().coord(i),
                        histogram.second().coord(j),
                        histogram(i, j)
                    );
                }
            }
        }

        if (options.autocorrelation) {
            for (auto& bond: bonds) {
                auto it = existing_bonds.find(bond);
                if (it == existing_bonds.end()) {
                    // New bond. Insert it and pad with zeros
                    auto pair = existing_bonds.emplace(bond, std::vector<float>(used_steps, 0.0));
                    pair.first->second.push_back(1.0);
                } else {
                    // Already seen this bond, add a single 1
                    it->second.push_back(1.0);
                }
            }
            // Add 0 to all bonds we did not see in this frame
            for (auto& it: existing_bonds) {
                if (it.second.size() != used_steps + 1) {
                    it.second.push_back(0.0);
                }
            }
        }
        used_steps += 1;
    }

    if (options.autocorrelation && used_steps != 0) {
        // Compute the autocorrelation for all bonds and average them
        auto correlator = Autocorrelation(used_steps);
        for (auto&& it: std::move(existing_bonds)) {
            correlator.add_timeserie(std::move(it.second));
        }
        correlator.normalize();
        auto& correlation = correlator.get_result();

        std::ofstream outcorr(options.autocorr_output, std::ios::out);
        if (!outcorr.is_open()) {
            throw CFilesError("Could not open the '" + options.autocorr_output + "' file.");
        }
        fmt::print(outcorr, "# Auto correlation between H-bonds existence\n");
        fmt::print(outcorr, "# step value\n");

        auto norm = correlation[0];
        for (size_t i=0; i<correlation.size() / 2; i++) {
            fmt::print(outcorr, "{} {}\n", i * options.steps.stride(), correlation[i] / norm);
        }
    }

    return 0;
}
