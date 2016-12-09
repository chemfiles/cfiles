// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#include <iostream>
#include <iterator>

#include "chemfiles.hpp"

#include "CommandFactory.hpp"
#include "commands/Help.hpp"

static auto CFILES_VERSION = "0.0.1";

static std::string version() {
    return std::string("version ") + CFILES_VERSION + " (using chemfiles " + CHEMFILES_VERSION + ")";
}

static void print_usage() {
    std::cout << "Usage: cfiles subcommand [--options] [args]" << std::endl;
    std::cout << std::endl;

    auto help = Help();
    help.list_commands();
}

int main(int argc, const char* argv[]) {
    if (argc < 2){
        std::cout << "cfiles: analysis algorithms for theoretical chemistry" << std::endl;
        std::cout << version() << std::endl;
        std::cout << std::endl;
        print_usage();
        return EXIT_FAILURE;
    }

    // Check first for version or help flags
    if (argv[1] == std::string("-h") || argv[1] == std::string("--help")) {
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
        std::cout << "cfiles " << version() << std::endl;
        return EXIT_SUCCESS;
    }

    // Error on other flags
    if (argv[1][0] == '-') {
        std::cout << "Unkown flag: " << argv[1] << std::endl;
        std::cout << std::endl;
        print_usage();
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
