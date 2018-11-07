// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_ROTATION_CORRELATION_HPP
#define CFILES_ROTATION_CORRELATION_HPP

#include <chemfiles.hpp>

#include "Command.hpp"
#include "utils.hpp"

class Rotcf final: public Command {
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
        /// Output file path
        std::string outfile;
        /// Selection for the orientation vector
        std::string selection;
    };

    Rotcf() {}
    int run(int argc, const char* argv[]) override;
    std::string description() const override;
};

#endif
