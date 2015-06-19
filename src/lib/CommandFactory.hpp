/*
 * chrp, an analysis frontend for the Chemharp library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#pragma once
#ifndef CHRP_FRONTEND_COMMAND_FACTORY_HPP
#define CHRP_FRONTEND_COMMAND_FACTORY_HPP

#include <memory>
#include <string>
#include <map>

#include "Command.hpp"
#include "Errors.hpp"

#include "Help.hpp"
#include "Rdf.hpp"

using command_creator_t = std::unique_ptr<Command>(*)(void);

static std::map<std::string, command_creator_t> COMMANDS = {
    {"help", [](){return std::unique_ptr<Command>(new Help());}},
    {"rdf", [](){return std::unique_ptr<Command>(new Rdf());}},
};

inline std::unique_ptr<Command> get_command(const std::string& name) {
    auto it = COMMANDS.find(name);
    if (it == COMMANDS.end()){
        throw chrp_exception("Can not find the subcommand '" + name + "'");
    }
    return it->second();
}

#endif
