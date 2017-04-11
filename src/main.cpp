// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <iostream>

#include "CommandFactory.hpp"
#include "utils.hpp"

static void list_commands();
static void print_usage();

int main(int argc, const char* argv[]) {
    // Check first for version or help flags
    if (argc < 2 || argv[1] == std::string("-h") || argv[1] == std::string("--help")) {
        print_usage();
        return 0;
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


static void print_usage() {
    std::cout << "cfiles: file algorithms for theoretical chemistry\n";
    std::cout << full_version() << std::endl;
    std::cout << R"(Guillaume Fraux <guillaume@fraux.fr>

cfiles provides file handling and analysis algorithms for theoretical chemistry
trajectories. Each algorithm is accessible through a specific sub-command: for
example `cfiles merge` to merge files; `cfiles rdf` for radial distribution
functions; etc. Input, output and selection capacities are provided by the
chemfiles library (http://chemfiles.org).

Use 'cfiles <command> --help' to get more information about a specific command.

Usage:
  cfiles <command> [--options] [args]

Examples:
  cfiles merge --help
  cfiles rdf water.tng -s "name O" --max=8.5 --output=rdf-O-O.dat
  cfiles angles result.xtc --topology=initial.mol --topology-format=PDB

)";

    list_commands();
}

static void list_commands() {
    const size_t command_width = 10;
    std::cout << "Available commands:" << std::endl;
    for (auto& command: all_commands()) {
        auto name = command.name;
        auto description = command.create()->description();
        std::cout << "  " <<  name << std::string(command_width - name.size(), ' ');
        std::cout << description << std::endl;
    }
}
