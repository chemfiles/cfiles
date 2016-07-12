/* cfiles, an analysis frontend for the Chemfiles library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#ifndef CFILES_AVERAGE_COMMAND_HPP
#define CFILES_AVERAGE_COMMAND_HPP

#include <map>
#include <chemfiles.hpp>

#include "Histogram.hpp"
#include "Command.hpp"

namespace docopt {
    struct value;
}

/// Base class for time-averaged computations
class AverageCommand: public Command {
public:
    struct Options {
        //! Input trajectory
        std::string trajectory;
        //! First step to use
        size_t start;
        //! Last step to use
        size_t end;
        //! Use a step every `stride` steps
        size_t stride;
        //! Do we have a custom cell to use?
        bool custom_cell;
        //! Unit cell to use
        chemfiles::UnitCell cell;
        //! Topology file to use
        std::string topology;
        //! Should we try to guess the topology?
        bool guess_bonds;
    };

    //! A strinc containing Doctopt style options for all time-averaged commands.
    //! It should be added to the command-specific options.
    static const std::string AVERAGE_OPTIONS;

    virtual ~AverageCommand() = default;
    int run(int argc, const char* argv[]) override final;

    //! Setup the command. This function should call `AverageCommand::parse_options`
    virtual void setup(int argc, const char* argv[], Histogram<double>& histogram) = 0;
    //! Add the data from a `frame` to the `histogram`
    virtual void accumulate(const chemfiles::Frame& frame, Histogram<double>& histogram) = 0;
    //! Finish the run, and write any output
    virtual void finish(const Histogram<double>& histogram) = 0;

protected:
    //! Get access to the options for this run
    const Options& options() const {return options_;}
    //! Parse the options from a doctop map/
    void parse_options(const std::map<std::string, docopt::value>& args);

private:
    //! Options for this instance of RDF
    Options options_;
    //! Histogram, for inserting the data at each step
    Histogram<double> histogram_;
    //! Result for storing the pre-normalized results
    std::vector<double> result_;
    //! Number of steps for the average
    size_t nsteps_ = 0;
};

#endif
