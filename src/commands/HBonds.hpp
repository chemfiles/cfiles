// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_HBONDS_HPP
#define CFILES_HBONDS_HPP

#include <chemfiles.hpp>

#include "Command.hpp"
#include "utils.hpp"

class HBonds final: public Command {
public:
    struct Options {
        /// Input trajectory
        std::string trajectory;
        /// Specific format to use with the trajectory
        std::string format;
        /// Specific steps to use from the trajectory
        steps_range steps;
        /// Do we have a custom cell to use?
        bool custom_cell = false;
        /// Unit cell to use
        chemfiles::UnitCell cell;
        /// Topology file to use
        std::string topology;
        /// Format to use for the topology file
        std::string topology_format;
        /// Should we try to guess the topology?
        bool guess_bonds = false;
        /// HBonds output
        std::string outfile;
        /// Should we compute the autocorrelation
        bool autocorrelation = false;
        /// Autocorrelation output
        std::string autocorr_output;
        /// Should we compute the hydrogen bonds histogram
        bool histogram = false;
        /// Autocorrelation output
        std::string histogram_output;
        /// Selection for the acceptor of the hydrogen bond (usually O/N/S)
        std::string acceptor_selection;
        /// Selection for the donor of the hydrogen bond (usually O-H/N-H)
        std::string donor_selection;
        /// Parameters for donor-acceptor max distance (in angstroms)
        double distance;
        /// and for acceptor-donor-hydrogen max angle (in degrees)
        double angle;
        /// If computing the histogram, how many points should it have
        size_t npoints;
    };

    HBonds() {}
    int run(int argc, const char* argv[]) override;
    std::string description() const override;
};

#endif
