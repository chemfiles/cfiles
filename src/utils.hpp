/* cfiles, an analysis frontend for the Chemfiles library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#pragma once
#ifndef CFILES_UTILS_HPP
#define CFILES_UTILS_HPP

#include <vector>
#include <string>
#include <sstream>

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

inline std::vector<double> parse_cell(const std::string& cell_string) {
    auto cell = std::vector<double>();
    auto cell_strs = split(cell_string, ':');
    if (cell_strs.size()==1) {
        cell.push_back(stod(cell_strs[0]));
        cell.push_back(stod(cell_strs[0]));
        cell.push_back(stod(cell_strs[0]));
    } else if (cell_strs.size()==3) {
        cell.push_back(stod(cell_strs[0]));
        cell.push_back(stod(cell_strs[1]));
        cell.push_back(stod(cell_strs[2]));
    } else if (cell_strs.size()==6) {
        cell.push_back(stod(cell_strs[0]));
        cell.push_back(stod(cell_strs[1]));
        cell.push_back(stod(cell_strs[2]));
        cell.push_back(stod(cell_strs[3]));
        cell.push_back(stod(cell_strs[4]));
        cell.push_back(stod(cell_strs[5]));
    } else {
        throw CFilesError("The cell string should have 1, 3 or 6 values.");
    }

    return cell;
}

#endif
