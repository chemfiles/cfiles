// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <sstream>

#include <chemfiles.hpp>
#include <chemfiles.h>

#include "version.hpp"
#include "utils.hpp"
#include "Errors.hpp"

double string2double(const std::string& string) {
    try {
        size_t length = 0;
        double value = std::stod(string, &length);
        if (length != string.length()) {
            throw cfiles_error("Can not convert '{}' to a number", string);
        }
        return value;
    } catch (const std::invalid_argument&) {
        throw cfiles_error("Can not convert '{}' to a number", string);
    }
}

long long string2long(const std::string& string) {
    try {
        size_t length = 0;
        long long value = std::stoll(string, &length);
        if (length != string.length()) {
            throw cfiles_error("Can not convert '{}' to a number", string);
        }
        return value;
    } catch (const std::invalid_argument&) {
        throw cfiles_error("Can not convert '{}' to a number", string);
    }
}

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

std::string trim(const std::string& str) {
    auto front = std::find_if_not(str.begin(), str.end(), [](int c) {
        return std::isspace(c);
    });
    auto back = std::find_if_not(str.rbegin(), str.rend(), [](int c) {
        return std::isspace(c);
    }).base();
    return (back <= front ? std::string() : std::string(front, back));
}


chemfiles::UnitCell parse_cell(const std::string& string) {
    auto splitted = split(string, ':');
    if (splitted.size() == 1) {
        auto a = string2double(splitted[0]);
        if (a <= 0) {
            throw cfiles_error("custom cell can not have negative length");
        }
        return chemfiles::UnitCell(a);
    } else if (splitted.size() == 3) {
        auto a = string2double(splitted[0]);
        auto b = string2double(splitted[1]);
        auto c = string2double(splitted[2]);
        if (a <= 0 || b <= 0 || c <= 0) {
            throw cfiles_error("custom cell can not have negative length");
        }
        return chemfiles::UnitCell(a, b, c);
    } else if (splitted.size() == 6) {
        auto a = string2double(splitted[0]);
        auto b = string2double(splitted[1]);
        auto c = string2double(splitted[2]);
        auto alpha = string2double(splitted[3]);
        auto beta  = string2double(splitted[4]);
        auto gamma = string2double(splitted[5]);
        if (a <= 0 || b <= 0 || c <= 0) {
            throw cfiles_error("custom cell can not have negative length");
        }
        if (alpha <= 0 || beta <= 0 || gamma <= 0) {
            throw cfiles_error("custom cell can not have negative angles");
        }
        return chemfiles::UnitCell(a, b, c, alpha, beta, gamma);
    } else {
        throw cfiles_error(
            "custom cell should be specified as 'a:b:c:α:β:γ'' or 'a:b:c' or 'a'"
        );
    }
}

steps_range steps_range::parse(const std::string& string) {
    steps_range range;
    auto splitted = split(string, ':');
    if (splitted.size() == 3 || splitted.size() == 2) {
        if (splitted[0] != "") {
            auto first = string2long(splitted[0]);
            if (first < 0) {
                throw cfiles_error(
                    "starting step must be positive, not {}", first
                );
            } else {
                range.first_ = first;
            }
        }

        if (splitted[1] != "") {
            auto last = string2long(splitted[1]);
            if (last < 0) {
                throw cfiles_error(
                    "last step must be positive, not {}", last
                );
            } else {
                range.last_ = last;
            }
        }

        if (range.last_ < range.first_) {
            throw cfiles_error(
                "last step ({}) must be bigger than the first step ({})",
                range.last_, range.first_
            );
        }

        if (splitted.size() == 3) {
            if (splitted[2] != "") {
                auto stride = string2long(splitted[2]);
                if (stride <= 0) {
                    throw cfiles_error(
                        "stride must be bigger than 1, not {}", stride
                    );
                } else {
                    range.stride_ = stride;
                }
            }
        }
    } else {
        throw cfiles_error(
            "steps range should be specified as 'start:stop' or 'start:stop:stride'"
        );
    }

    if (range.last_ != static_cast<size_t>(-1)) {
        if (range.last_ % range.stride_ != 0) {
            auto div = range.last_ / range.stride_;
            range.last_ = (div + 1) * range.stride_;
        }
    }
    return range;
}
