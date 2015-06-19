/*
 * chrp, an analysis frontend for the Chemharp library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#pragma once
#ifndef CHRP_FRONTEND_COMMAND_HPP
#define CHRP_FRONTEND_COMMAND_HPP

#include <string>

/**
 * @class Command Command.hpp Command.cpp
 *
 * Basic subcommand for chrp. The only method is run, which will be called with the
 * arguments of the subcommand. If the command is:
 * 		chrp sub --here -i -k -lkj positional positional -o positional
 * Then, the run function will be called with
 * 	    args = ["sub", "--here", "-i", "-k", "-lkj", "positional", "positional", "-o", "positional"]
 */
class Command {
public:
    Command() = default;
    virtual ~Command() = default;
    //! Run the command, with the arguments \c argv, of size argc.
    virtual int run(int argc, char** argv) = 0;
    //! Basic description of the command
    virtual std::string description() const = 0;
    //! More detailed help of the command
    virtual std::string help() const {
        return description();
    }
};

#endif
