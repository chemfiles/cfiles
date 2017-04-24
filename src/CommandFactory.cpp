// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include "CommandFactory.hpp"

#include "commands/Angles.hpp"
#include "commands/Convert.hpp"
#include "commands/Merge.hpp"
#include "commands/Rdf.hpp"
#include "commands/HBonds.hpp"
#include "commands/DensityProfile.hpp"

const std::vector<command_creator>& all_commands() {
    static std::vector<command_creator> commands = {
        {"angles", [](){return std::unique_ptr<Command>(new AngleDistribution());}},
        {"convert", [](){return std::unique_ptr<Command>(new Convert());}},
        {"merge", [](){return std::unique_ptr<Command>(new Merge());}},
        {"rdf", [](){return std::unique_ptr<Command>(new Rdf());}},
        {"hbonds", [](){return std::unique_ptr<Command>(new HBonds());}},
        {"density", [](){return std::unique_ptr<Command>(new DensityProfile());}},
    };
    return commands;
}
