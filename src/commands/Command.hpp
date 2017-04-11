// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_COMMAND_HPP
#define CFILES_COMMAND_HPP

#include <string>


/// Basic subcommand for `cfiles`. The main method is `run`, which will be called
/// with the arguments of the subcommand.
class Command {
public:
    Command() = default;
    virtual ~Command() = default;

    /// Run the command, with the arguments array `argv`, containing `argc`
    /// strings. This should output help to the standard output when `argv`
    /// contains the string `"--help"`.
    virtual int run(int argc, const char* argv[]) = 0;
    /// Output a description of the command
    virtual std::string description() const = 0;
};

#endif
