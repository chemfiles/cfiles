// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_ELASTIC_HPP
#define CFILES_ELASTIC_HPP

#include "Command.hpp"
#include "utils.hpp"

namespace chemfiles {
    class Matrix3D;
}

class Elastic final: public Command {
public:
    struct Options {
        /// Input trajectory
        std::string trajectory;
        /// Specific format to use with the trajectory
        std::string format = "";
        /// Specific steps to use from the trajectory
        steps_range steps;
        /// Output data file
        std::string outfile;
        /// Temperature of the simulation
        double temperature;
    };

    Elastic() {}
    std::string description() const override;
    int run(int argc, const char* argv[]) override;

private:
};

#endif
