// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_COMMAND_FACTORY_HPP
#define CFILES_COMMAND_FACTORY_HPP

#include <memory>
#include <vector>

#include "commands/Command.hpp"
#include "Errors.hpp"

struct command_creator {
    using command_creator_t = std::unique_ptr<Command>(*)(void);
    /// Command name
    std::string name;
    /// Command instanciation
    command_creator_t create;
};

const std::vector<command_creator>& all_commands();

inline std::unique_ptr<Command> get_command(const std::string& name) {
    for (auto& command: all_commands()) {
        if (command.name == name) {
            return command.create();
        }
    }
    throw cfiles_error("No subcommand named '{}' available", name);
}

#endif
