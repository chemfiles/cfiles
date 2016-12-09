// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#ifndef CFILES_UTILS_HPP
#define CFILES_UTILS_HPP

#include <chemfiles.hpp>

#include "Errors.hpp"

inline std::vector<std::string> split(const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> elems;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

inline chemfiles::UnitCell parse_cell(const std::string& cell_string) {
    auto splitted = split(cell_string, ':');
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

#endif
