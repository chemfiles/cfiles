// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

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
