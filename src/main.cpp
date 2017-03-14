// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#include <iostream>

#include "CommandFactory.hpp"
#include "commands/Help.hpp"
#include "utils.hpp"

int main(int argc, const char* argv[]) {
    // Check first for version or help flags
    if (argc < 2 || argv[1] == std::string("-h") || argv[1] == std::string("--help")) {
        auto command = get_command("help");
        const char* args[] = {""};
        return command->run(1, args);
    }

    auto find_in_argv = [&argc, &argv](const std::string& value) -> bool {
        for (int i = 0; i<argc; i++) {
            if (argv[i] == value) {
                return true;
            }
        }
        return false;
    };

    if (find_in_argv("-V") || find_in_argv("--version")) {
        std::cout << "cfiles " << full_version() << std::endl;
        return EXIT_SUCCESS;
    }

    // Error on other flags
    if (argv[1][0] == '-') {
        std::cout << "Unkown flag: " << argv[1] << std::endl;
        return EXIT_FAILURE;
    }

    try {
        auto command = get_command(argv[1]);
        return command->run(argc - 1, &argv[1]);
    } catch (const std::exception& e){
        std::cout << "Error: " << e.what() << std::endl;
        return 2;
    }
}
