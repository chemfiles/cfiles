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
        return std::stod(string);
    } catch (const std::invalid_argument&) {
        throw CFilesError("Can not convert '" + string + "' to number");
    }
}

long string2long(const std::string& string) {
    try {
        return std::stol(string);
    } catch (const std::invalid_argument&) {
        throw CFilesError("Can not convert '" + string + "' to number");
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

chemfiles::UnitCell parse_cell(const std::string& string) {
    auto splitted = split(string, ':');
    if (splitted.size() == 1) {
        return chemfiles::UnitCell(string2double(splitted[0]));
    } else if (splitted.size() == 3) {
        auto a = string2double(splitted[0]);
        auto b = string2double(splitted[1]);
        auto c = string2double(splitted[2]);
        return chemfiles::UnitCell(a, b, c);
    } else if (splitted.size() == 6) {
        auto a = string2double(splitted[0]);
        auto b = string2double(splitted[1]);
        auto c = string2double(splitted[2]);
        auto alpha = string2double(splitted[3]);
        auto beta  = string2double(splitted[4]);
        auto gamma = string2double(splitted[5]);
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
            long long first = -1;
            try {
                first = std::stoll(splitted[0]);
            } catch(const std::invalid_argument&) {
                throw CFilesError("Can not convert '" + string + "' to number");
            }
            if (first < 0) {
                throw CFilesError(
                    "starting step must be positive, not " + std::to_string(first)
                );
            } else {
                range.first_ = first;
            }
        }

        if (splitted[1] != "") {
            long long last = -1;
            try {
                last = std::stoll(splitted[1]);
            } catch(const std::invalid_argument&) {
                throw CFilesError("Can not convert '" + string + "' to number");
            }
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
                long long stride = -1;
                try {
                    stride = std::stoll(splitted[2]);
                } catch(const std::invalid_argument&) {
                    throw CFilesError("Can not convert '" + string + "' to number");
                }
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
        // integer division to know the number of steps to use
        auto div = range.last_ / range.stride_;
        // last step to use (useful if last % stride !=0)
        range.last_ = div * range.stride_;
    }
    return range;
}
