// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_COMMAND_HPP
#define CFILES_COMMAND_HPP

#include <Options.hpp>

/// Basic subcommand for `cfiles`.
class Command {
public:
    Command(Options options): options_(std::move(options)) {}
    virtual ~Command() = default;

    /// Get the command options, the user should call `Options::parse` with the
    /// command line argument, or set any required option manually.
    Options& options() {
        return options_;
    }

    /// Run the command, the options should have been set before calling this.
    virtual int run() = 0;

protected:
    Options options_;
};

#endif
