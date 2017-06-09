// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <iostream>
#include <set>
#include "warnings.hpp"

void warn(std::string message) {
    std::cerr << "[cfiles] " << message << std::endl;
}

void warn_once(std::string message) {
    static std::set<std::string> ALREADY_SEEN;
    auto not_seen = ALREADY_SEEN.insert(message).second;
    if (not_seen) {
        warn(message);
    }
}
