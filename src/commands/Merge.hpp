// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#ifndef CFILES_MERGE_HPP
#define CFILES_MERGE_HPP

#include <chemfiles.hpp>

#include "Command.hpp"

class Merge final: public Command {
public:
    struct Options {
        std::vector<std::string> infiles;
        std::vector<std::string> input_formats;
        std::string outfile;
        std::string output_format = "";
        bool custom_cell = false;
        chemfiles::UnitCell cell;
    };

    Merge() {}
    int run(int argc, const char* argv[]) override;
    std::string description() const override;

private:
    Options options_;
};

#endif
