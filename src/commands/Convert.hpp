// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#ifndef CFILES_CONVERT_HPP
#define CFILES_CONVERT_HPP

#include <chemfiles.hpp>

#include "Command.hpp"

class Convert final: public Command {
public:
    struct Options {
        std::string infile;
        std::string outfile;
        std::string input_format = "";
        std::string output_format = "";
        std::string topology = "";
        std::string topology_format = "";
        bool custom_cell = false;
        chemfiles::UnitCell cell;
        bool guess_bonds = false;
        bool wrap = false;
    };

    Convert() {}
    int run(int argc, const char* argv[]) override;
    std::string description() const override;
private:
    Options options_;
};

#endif
