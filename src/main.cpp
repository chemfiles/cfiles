/* cfiles, an analysis frontend for the Chemfiles library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#include <iostream>

#include "CommandFactory.hpp"
#include "Help.hpp"


int usage() {
    std::cout << "Usage: cfiles subcommand [--options] [args]" << std::endl;
    std::cout << std::endl;

    auto help = Help();
    help.list_commands();

    return 1;
}

int main(int argc, const char* argv[]) {
    if (argc < 2){
        return usage();
    }

    try {
        auto command = get_command(argv[1]);
        return command->run(argc - 1, &argv[1]);
    } catch (const std::exception& e){
        std::cout << "Error: " << e.what() << std::endl;
        return 2;
    }
}
