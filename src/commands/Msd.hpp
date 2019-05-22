// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_MSD_HPP
#define CFILES_MSD_HPP

#include <chemfiles.hpp>

#include "Command.hpp"
#include "utils.hpp"

class MSD final: public Command {
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
        /// msd output
        std::string outfile;
        /// Selection of atoms to use when computing MSD
        std::string selection;
        /// Should we unwrap the positions?
        bool unwrap = false;
    };

    MSD() {}
    int run(int argc, const char* argv[]) override;
    std::string description() const override;
};

#endif
