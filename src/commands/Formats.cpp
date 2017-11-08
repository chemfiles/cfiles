// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <iostream>

#include <fmt/format.h>
#include <docopt/docopt.h>
#include <chemfiles/FormatFactory.hpp>

#include "Formats.hpp"
#include "utils.hpp"

using namespace chemfiles;

static const char OPTIONS[] =
R"(List available file formats

Usage:
  cfiles formats
  cfiles formats (-h | --help)

Examples:
    cfiles formats

Options:
  -h --help                     show this help
)";

static void parse_options(int argc, const char* argv[]) {
    auto options_str = command_header("info", Formats().description());
    options_str += "Guillaume Fraux <guillaume@fraux.fr>\n\n";
    options_str += OPTIONS;
    auto args = docopt::docopt(options_str, {argv, argv + argc}, true, "");
}

std::string Formats::description() const {
    return "list available file formats";
}

int Formats::run(int argc, const char* argv[]) {
    parse_options(argc, argv);

    fmt::print("Available formats [name (extension) description]:\n\n");
    for (auto format: FormatFactory::get().formats()) {
        auto name = format.name();
        auto ext = format.extension();
        auto description = format.description();

        int fill = 18 - name.length() - ext.length();
        if (fill < 0) { fill = 0; }
        fmt::print("{} ({}) {:{}}{} \n", name, ext, " ", fill, description);
    }

    return 0;
}
