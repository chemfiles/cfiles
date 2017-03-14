// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#ifndef CFILES_UTILS_HPP
#define CFILES_UTILS_HPP

#include <string>
#include <vector>

namespace chemfiles {
    class UnitCell;
}

/// Get a string describing the full version of cfiles
std::string full_version();

/// Create the command description header
std::string command_header(std::string name, std::string description);

/// Split a string a delimiter
std::vector<std::string> split(const std::string& string, char delimiter);

/// Parse an unit cell string
chemfiles::UnitCell parse_cell(const std::string& string);

#endif
