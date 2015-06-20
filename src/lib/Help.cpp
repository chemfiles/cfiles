/*
 * chrp, an analysis frontend for the Chemharp library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#include "Help.hpp"
#include "CommandFactory.hpp"

std::string Help::description() const{
    return "Get help about the subcommands.";
}

int Help::run(int argc, char** argv){
    if (argc == 1) {
        list_commands();
    } else {
        about(argv[1]);
    }
    return 0;
}

void Help::list_commands() const {
    std::cout << "Use 'chrp help <subcommand>' to get help about a specific subcommand." ;
    std::cout << std::endl << std::endl;

    const size_t command_width = 10;
    std::cout << "Available subcommands:" << std::endl;
    for (auto it : COMMANDS()) {
        auto name = it.first;
        auto command = it.second();
        std::cout << "  " << name  << std::string(command_width - name.size(), ' ');
        std::cout << command->description() << std::endl;
    }
}

void Help::about(const std::string& name) const {
    auto command = get_command(name);
    std::cout << command->help() << std::endl;
}
