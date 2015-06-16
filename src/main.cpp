/*
 * chrp, an analysis frontend for the Chemharp library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#include <iostream>
#include <map>
#include <memory>

#include "lib/Command.hpp"

int usage(const char* name);

using command_creator_t = std::unique_ptr<Command>(*)(void);
static std::map<std::string, command_creator_t> COMMANDS = {
};


int main(int argc, char** argv) {
    if (argc < 2){
        return usage(argv[0]);
    }

    auto com = COMMANDS.find(argv[1]);
    int res = 0;
    if (com != COMMANDS.end()){
        auto command = (*com).second();
        try {
            res = command->run(argc - 1, &argv[1]);
        } catch (const std::exception& e){
            std::cout << "Error: " << e.what() << std::endl;
            res = -1;
        }
    } else {
        std::cout << "Could not find the '" << argv[1] << "' subcommand." << std::endl;
        std::cout << std::endl;
        res = usage(argv[0]);
    }
    return res;
}

int usage(const char* name) {
    std::cout << "usage: " << name << " subcommand [--options] args" << std::endl;
    return 1;
}
