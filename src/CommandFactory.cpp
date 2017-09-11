// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include "CommandFactory.hpp"

#include "commands/Angles.hpp"
#include "commands/Convert.hpp"
#include "commands/DensityProfile.hpp"
#include "commands/Merge.hpp"
#include "commands/HBonds.hpp"
#include "commands/Info.hpp"
#include "commands/Rdf.hpp"

const std::vector<command_creator>& all_commands() {
    static std::vector<command_creator> commands = {
        {"angles", [](){return std::unique_ptr<Command>(new AngleDistribution());}},
        {"convert", [](){return std::unique_ptr<Command>(new Convert());}},
        {"density", [](){return std::unique_ptr<Command>(new DensityProfile());}},
        {"hbonds", [](){return std::unique_ptr<Command>(new HBonds());}},
        {"info", [](){return std::unique_ptr<Command>(new Info());}},
        {"merge", [](){return std::unique_ptr<Command>(new Merge());}},
        {"rdf", [](){return std::unique_ptr<Command>(new Rdf());}},
    };
    return commands;
}
