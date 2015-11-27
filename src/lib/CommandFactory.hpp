/* cfiles, an analysis frontend for the Chemfiles library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#pragma once
#ifndef CFILES_COMMAND_FACTORY_HPP
#define CFILES_COMMAND_FACTORY_HPP

#include <memory>
#include <string>
#include <map>

#include "Command.hpp"
#include "Errors.hpp"

using command_creator_t = std::unique_ptr<Command>(*)(void);
const std::map<std::string, command_creator_t>& COMMANDS();

inline std::unique_ptr<Command> get_command(const std::string& name) {
    auto it = COMMANDS().find(name);
    if (it == COMMANDS().end()){
        throw CFilesError("Can not find the subcommand '" + name + "'");
    }
    return it->second();
}

#endif
