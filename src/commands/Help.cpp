// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#include <iostream>
#include <chemfiles.h>

#include "Help.hpp"
#include "CommandFactory.hpp"
#include "utils.hpp"

std::string Help::description() const{
    return "get help about subcommands";
}

int Help::run(int argc, const char* argv[]){
    if (argc == 1) {
        std::cout << "cfiles: analysis algorithms for theoretical chemistry" << std::endl;
        std::cout << full_version() << std::endl;
        std::cout << "Guillaume Fraux <guillaume@fraux.fr>" << std::endl << std::endl;

        std::cout << "Usage:" << std::endl;
        std::cout << "  cfiles <command> [--options] [args]" << std::endl << std::endl;

        std::cout << "Use 'cfiles help <command>' to get help about a specific command." << std::endl << std::endl;
        list_commands();
    } else {
        about(argv[1]);
    }
    return 0;
}

void Help::list_commands() const {
    const size_t command_width = 10;
    std::cout << "Available commands:" << std::endl;
    for (auto& command: all_commands()) {
        auto name = command.name;
        auto description = command.create()->description();
        std::cout << "  " <<  name << std::string(command_width - name.size(), ' ');
        std::cout << description << std::endl;
    }
}

void Help::about(const std::string& name) const {
    auto command = get_command(name);
    const char* argv[3] = {"cfiles", name.c_str(), "--help"};
    std::cout << command->run(3, argv) << std::endl;
}
