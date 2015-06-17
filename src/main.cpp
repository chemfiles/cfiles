/*
 * chrp, an analysis frontend for the Chemharp library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#include <iostream>

#include "lib/CommandFactory.hpp"

int usage() {
    std::cout << "Usage: chrp subcommand [--options] [args]" << std::endl;
    std::cout << std::endl;
    std::cout << "Some usefull subcommands:" << std::endl;

    const std::string SEP = "    ";
    auto help = get_command("help");
    std::cout << SEP << "'help' " << help->description() << std::endl;

    return 1;
}

int main(int argc, char** argv) {
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
