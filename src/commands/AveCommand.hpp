// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#ifndef CFILES_AVERAGE_COMMAND_HPP
#define CFILES_AVERAGE_COMMAND_HPP

#include <map>
#include <chemfiles.hpp>

#include "Averager.hpp"
#include "Command.hpp"
#include "utils.hpp"

namespace docopt {
    struct value;
}

/// Base class for time-averaged computations
class AveCommand: public Command {
public:
    struct Options {
        /// Input trajectory
        std::string trajectory;
        /// Specific format to use with the trajectory
        std::string format = "";
        /// Specific steps to use from the trajectory
        steps_range steps;
        /// Do we have a custom cell to use?
        bool custom_cell;
        /// Unit cell to use
        chemfiles::UnitCell cell;
        /// Topology file to use
        std::string topology = "";
        /// Format to use for the topology file
        std::string topology_format = "";
        /// Should we try to guess the topology?
        bool guess_bonds = false;
    };

    /// A strinc containing Doctopt style options for all time-averaged commands.
    /// It should be added to the command-specific options.
    static const std::string AVERAGE_OPTIONS;

    virtual ~AveCommand() = default;
    int run(int argc, const char* argv[]) override final;

    /// Setup the command and the histogram.
    /// This function MUST call `AverageCommand::parse_options`.
    virtual Averager<double> setup(int argc, const char* argv[]) = 0;
    /// Add the data from a `frame` to the `histogram`
    virtual void accumulate(const chemfiles::Frame& frame, Histogram<double>& histogram) = 0;
    /// Finish the run, and write any output
    virtual void finish(const Histogram<double>& histogram) = 0;

protected:
    /// Get access to the options for this run
    const Options& options() const {return options_;}
    /// Parse the options from a doctop map/
    void parse_options(const std::map<std::string, docopt::value>& args);

private:
    /// Options
    Options options_;
    /// Averaging histogram for the data
    Averager<double> histogram_;
};

#endif
