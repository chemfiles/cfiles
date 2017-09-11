// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_INFO_HPP
#define CFILES_INFO_HPP

#include "Command.hpp"

class Info final: public Command {
public:
    struct Options {
        std::string input;
        bool guess_bonds;
        size_t step;
    };

    Info() {}
    int run(int argc, const char* argv[]) override;
    std::string description() const override;
};

#endif
