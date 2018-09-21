// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <fstream>

#include <docopt/docopt.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <chemfiles.hpp>
#include <Eigen/Dense>

#include "Elastic.hpp"
#include "Errors.hpp"

using namespace chemfiles;

using Matrix6 = Eigen::Matrix<double, 6, 6>;

static size_t CARTESIAN_TO_VOIGT[][2] = {
    {0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1},
};

static const std::string OPTIONS =
R"(Compute the elastic tensor of a system from the unit cell fluctuations during
a NPT simulation.

The values given here are highly dependent on having good statistic during the
simulation: both in term of having a long enough simulation time to get to the
equilibrium, and using a good barostat that does produce isobaric-isothermal
ensemble fluctuations (not just average). The theory behind this code is
described in https://dx.doi.org/10.1080/08927022.2017.1313418.

Usage:
  cfiles elastic [options] <trajectory>
  cfiles elastic (-h | --help)

Examples:
  cfiles elastic -t 328 -o elastic.dat trajectory.pdb

Options:
  -h --help                        show this help
  --format=<format>                force the input file format to be <format>
  -t <temp>, --temperature=<temp>  temperature of the simulation, in kelvin
  -o <file>, --output=<file>       write result to <file>. This default to the
                                   trajectory file name with the `.angles.dat`
                                   extension.
  --steps=<steps>                  steps to use from the input. <steps> format
                                   is <start>:<end>[:<stride>] with <start>,
                                   <end> and <stride> optional. The used steps
                                   goes from <start> to <end> (excluded) by
                                   steps of <stride>. The default values are 0
                                   for <start>, the number of steps for <end>
                                   and 1 for <stride>.
)";


static Elastic::Options parse_options(int argc, const char* argv[]) {
    auto options_str = command_header("elastic", Elastic().description());
    options_str += "Guillaume Fraux <guillaume@fraux.fr>\n\n";
    options_str += std::string(OPTIONS);
    auto args = docopt::docopt(options_str, {argv, argv + argc}, true, "");

    Elastic::Options options;
    options.trajectory = args.at("<trajectory>").asString();

    if (args["--temperature"]) {
        options.temperature = string2double(args.at("--temperature").asString());
    } else {
        throw CFilesError("missing --temperature argument");
    }

    if (args.at("--steps")) {
        options.steps = steps_range::parse(args.at("--steps").asString());
    }

    if (args["--output"]) {
        options.outfile = args.at("--output").asString();
    } else {
        options.outfile = options.trajectory + ".elastic.dat";
    }

    return options;
}

std::string Elastic::description() const {
    return "compute elastic constants from unit cell fluctuations in NPT";
}

int Elastic::run(int argc, const char* argv[]) {
    auto options = parse_options(argc, argv);
    auto cells = std::vector<Matrix3D>();
    cells.reserve(options.steps.count());

    auto trajectory = Trajectory(options.trajectory, 'r', options.format);
    for (auto step: options.steps) {
        if (step >= trajectory.nsteps()) {
            break;
        }
        auto frame = trajectory.read_step(step);
        cells.emplace_back(frame.cell().matrix());
    }

    // Use the avrerage as the reference state
    auto reference = Matrix3D::zero();
    for (auto& cell: cells) {
        reference += cell.invert();
    }
    reference /= cells.size();
    auto reference_t = reference.transpose();

    auto epsilons = std::vector<Matrix3D>();
    epsilons.reserve(cells.size());
    for (auto& cell: cells) {
        epsilons.emplace_back(0.5 * (reference_t * cell.transpose() * cell * reference - Matrix3D::unit()));
    }

    // in GPa A^2 / K
    auto BOLTZMANN = 1.38065e-2;
    auto volume_inv = reference.determinant();
    auto v_kt = 1.0 / (volume_inv * BOLTZMANN * options.temperature);
    auto compute = [&](size_t ij[2], size_t kl[2]) {
        auto i = ij[0];
        auto j = ij[1];
        auto k = kl[0];
        auto l = kl[1];

        // Multiplicative factors for cross terms yz xz xy
        double factor = 1.0;
        if (i != j) {
            factor *= 2;
        }
        if (k != l) {
            factor *= 2;
        }

        double eij = 0;
        double ekl = 0;
        double eij_ekl = 0;
        for (auto& epsilon: epsilons) {
            eij += epsilon[i][j];
            ekl += epsilon[k][l];
            eij_ekl += epsilon[i][j] * epsilon[k][l];
        }
        eij /= epsilons.size();
        ekl /= epsilons.size();
        eij_ekl /= epsilons.size();

        return factor * v_kt * (eij_ekl - eij * ekl);
    };

    auto SVoigt = Matrix6();
    for (size_t i=0; i<6; i++) {
        for (size_t j=0; j<=i; j++) {
            SVoigt(i, j) = compute(CARTESIAN_TO_VOIGT[i], CARTESIAN_TO_VOIGT[j]);
        }
    }
    // Make the matrix symetric
    for (size_t i=0; i<6; i++) {
        for (size_t j=i+1; j<6; j++) {
            SVoigt(i, j) = SVoigt(j, i);
        }
    }

    if (std::abs(SVoigt.determinant()) < 100 * DBL_EPSILON) {
        throw CFilesError("the compliance matrix is not invertible");
    }

    auto CVoigt = SVoigt.inverse();
    std::ofstream outfile(options.outfile, std::ios::out);
    if (!outfile.is_open()) {
        throw CFilesError("Could not open the '" + options.outfile + "' file.");
    }

    fmt::print(outfile, "# stiffness tensor in GPa from {}\n", options.trajectory);
    for (size_t i=0; i<6; i++) {
        for (size_t j=0; j<6; j++) {
            if (j != 0) {
                fmt::print(outfile, " ");
            }
            if (j >= i) {
                fmt::print(outfile, "{:12.5f}", CVoigt(i, j));
            } else {
                fmt::print(outfile, "            ");
            }
        }
        fmt::print(outfile, "\n");
    }

    auto eigenvalues = CVoigt.eigenvalues();
    auto sorter = [](std::complex<double> i, std::complex<double> j) {
        return std::abs(i) < std::abs(j);
    };
    std::sort(eigenvalues.data(), eigenvalues.data() + eigenvalues.size(), sorter);
    fmt::print(outfile, "# eigen values of the stiffness tensor (GPa)\n");
    for (size_t i=0; i<6; i++) {
        auto& value = eigenvalues(i);
        if (value.imag() == 0) {
            fmt::print(outfile, "{:12.5f}\n", value.real());
        } else {
            fmt::print(outfile, "{:12.5f} + {:12.5f}i\n", value.real(), value.imag());
        }
    }

    double A = (CVoigt(0, 0) + CVoigt(1, 1) + CVoigt(2, 2)) / 3.0;
    double B = (CVoigt(1, 2) + CVoigt(0, 2) + CVoigt(0, 1)) / 3.0;
    double C = (CVoigt(3, 3) + CVoigt(4, 4) + CVoigt(5, 5)) / 3.0;

    double a = (SVoigt(0, 0) + SVoigt(1, 1) + SVoigt(2, 2)) / 3.0;
    double b = (SVoigt(1, 2) + SVoigt(0, 2) + SVoigt(0, 1)) / 3.0;
    double c = (SVoigt(3, 3) + SVoigt(4, 4) + SVoigt(5, 5)) / 3.0;

    double KV = (A + 2.0 * B) / 3.0;
    double GV = (A - B + 3.0 * C) / 5.0;
    double YV = 1.0 / (1.0 / (3.0 * GV) + 1.0 / (9.0 * KV));
    double PV = (1.0 - 3.0 * GV / (3.0 * KV + GV)) / 2.0;

    double KR = 1.0 / (3.0 * a + 6.0 * b);
    double GR = 5.0 / (4.0 * a - 4.0 * b + 3.0 * c);
    double YR = 1.0 / (1.0 / (3.0 * GR) + 1.0 / (9.0 * KR));
    double PR = (1.0 - 3.0 * GR / (3.0 * KR + GR)) / 2.0;

    double KH = (KV + KR) / 2.0;
    double GH = (GV + GR) / 2.0;
    double YH = 1.0 / (1.0 / (3.0 * GH) + 1.0 / (9.0 * KH));
    double PH = (1.0 - 3.0 * GH / (3.0 * KH + GH)) / 2.0;

    fmt::print(outfile, "# Bulk modulus (GPa) | Young's modulus (GPa) | Shear modulus (GPa) | Poisson's ratio\n");
    fmt::print(outfile, "# Voigt averaging\n");
    fmt::print(outfile, "{:12.5f} {:12.5f} {:12.5f} {:12.5f}\n", KV, YV, GV, PV);
    fmt::print(outfile, "# Reuss averaging\n");
    fmt::print(outfile, "{:12.5f} {:12.5f} {:12.5f} {:12.5f}\n", KR, YR, GR, PR);
    fmt::print(outfile, "# Hill averaging\n");
    fmt::print(outfile, "{:12.5f} {:12.5f} {:12.5f} {:12.5f}\n", KH, YH, GH, PH);

    return 0;
}
