// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#include "CommandFactory.hpp"

#include "commands/Angles.hpp"
#include "commands/Convert.hpp"
#include "commands/Merge.hpp"
#include "commands/Rdf.hpp"
#include "commands/HBonds.hpp"

const std::vector<command_creator>& all_commands() {
    static std::vector<command_creator> commands = {
        {"angles", [](){return std::unique_ptr<Command>(new AngleDistribution());}},
        {"convert", [](){return std::unique_ptr<Command>(new Convert());}},
        {"merge", [](){return std::unique_ptr<Command>(new Merge());}},
        {"rdf", [](){return std::unique_ptr<Command>(new Rdf());}},
        {"hbonds", [](){return std::unique_ptr<Command>(new HBonds());}},
    };
    return commands;
}
