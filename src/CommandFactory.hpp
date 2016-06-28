/* cfiles, an analysis frontend for the Chemfiles library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#ifndef CFILES_COMMAND_FACTORY_HPP
#define CFILES_COMMAND_FACTORY_HPP

#include <memory>
#include <string>
#include <vector>

#include "commands/Command.hpp"
#include "Errors.hpp"

struct command_creator {
    using command_creator_t = std::unique_ptr<Command>(*)(void);
    //! Command name
    std::string name;
    //! Command instanciation
    command_creator_t create;
};

const std::vector<command_creator>& all_commands();

inline std::unique_ptr<Command> get_command(const std::string& name) {
    for (auto& command: all_commands()) {
        if (command.name == name) {
            return command.create();
        }
    }
    throw CFilesError("Can not find the subcommand '" + name + "'");
}

#endif
