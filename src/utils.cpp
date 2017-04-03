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
    return header;
}

std::vector<std::string> split(const std::string& string, char delimiter) {
    std::stringstream sstream(string);
    std::string item;
    std::vector<std::string> tokens;
    while (std::getline(sstream, item, delimiter)) {
        tokens.push_back(item);
    }
    if (string.size() != 0 && string.back() == delimiter) {
        tokens.push_back("");
    }
    return tokens;
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
        throw CFilesError(
            "custom cell should be specified as 'a:b:c:α:β:γ'' or 'a:b:c' or 'a'"
        );
    }
}

steps_range steps_range::parse(const std::string& string) {
    steps_range range;
    auto splitted = split(string, ':');
    if (splitted.size() == 3 || splitted.size() == 2) {
        if (splitted[0] != "") {
            auto first = std::stoll(splitted[0]);
            if (first < 0) {
                throw CFilesError(
                    "starting step must be positive, not " + std::to_string(first)
                );
            } else {
                range.first_ = first;
            }
        }

        if (splitted[1] != "") {
            auto last = std::stoll(splitted[1]);
            if (last < 0) {
                throw CFilesError(
                    "last step must be positive, not " + std::to_string(last)
                );
            } else {
                range.last_ = last;
            }
        }

        if (range.last_ < range.first_) {
            throw CFilesError(
                "last step (" + std::to_string(range.last_) +
                ") must be bigger than the first step (" +
                std::to_string(range.first_) + ")"
            );
        }

        if (splitted.size() == 3) {
            if (splitted[2] != "") {
                auto stride = std::stoll(splitted[2]);
                if (stride <= 0) {
                    throw CFilesError(
                        "stride must be bigger than 1, not " + std::to_string(stride)
                    );
                } else {
                    range.stride_ = stride;
                }
            }
        }
    } else {
        throw CFilesError(
            "steps range should be specified as 'start:stop' or 'start:stop:stride'"
        );
    }

    if (range.last_ != static_cast<size_t>(-1)) {
        range.last_ = range.last_ / range.stride_ + range.stride_;
    }
    return range;
}
