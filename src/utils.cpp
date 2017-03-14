// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#include <sstream>

#include <chemfiles.hpp>
#include <chemfiles.h>

#include "version.hpp"
#include "utils.hpp"
#include "Errors.hpp"

std::string full_version() {
    return std::string("version ") + CFILES_VERSION + " (using chemfiles " + chfl_version() + ")";
}

std::string command_header(std::string name, std::string description) {
    auto header = std::string("cfiles ") + name + " (" + CFILES_VERSION + "): " + description + "\n";
    header += "Guillaume Fraux <guillaume@fraux.fr>\n";
    return header;
}

std::vector<std::string> split(const std::string& string, char delimiter) {
    std::stringstream ss(string);
    std::string item;
    std::vector<std::string> elems;
    while (std::getline(ss, item, delimiter)) {
        elems.push_back(item);
    }
    return elems;
}

chemfiles::UnitCell parse_cell(const std::string& string) {
    auto splitted = split(string, ':');
    if (splitted.size() == 1) {
        return chemfiles::UnitCell(stod(splitted[0]));
    } else if (splitted.size() == 3) {
        auto a = stod(splitted[0]);
        auto b = stod(splitted[1]);
        auto c = stod(splitted[2]);
        return chemfiles::UnitCell(a, b, c);
    } else if (splitted.size() == 6) {
        auto a = stod(splitted[0]);
        auto b = stod(splitted[1]);
        auto c = stod(splitted[2]);
        auto alpha = stod(splitted[3]);
        auto beta  = stod(splitted[4]);
        auto gamma = stod(splitted[5]);
        return chemfiles::UnitCell(a, b, c, alpha, beta, gamma);
    } else {
        throw CFilesError("The cell string should have 1, 3 or 6 values.");
    }
}
