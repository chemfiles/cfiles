// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_FORMATS_HPP
#define CFILES_FORMATS_HPP

#include "Command.hpp"

class Formats final: public Command {
public:
    Formats() {}
    int run(int argc, const char* argv[]) override;
    std::string description() const override;
};

#endif
