// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include "CommandFactory.hpp"

#include "commands/Angles.hpp"
#include "commands/Convert.hpp"
#include "commands/Density.hpp"
#include "commands/Elastic.hpp"
#include "commands/Formats.hpp"
#include "commands/HBonds.hpp"
#include "commands/Info.hpp"
#include "commands/Merge.hpp"
#include "commands/Msd.hpp"
#include "commands/Rdf.hpp"
#include "commands/Rotcf.hpp"

const std::vector<command_creator>& all_commands() {
    static std::vector<command_creator> commands = {
        {"angles", [](){return std::unique_ptr<Command>(new Angles());}},
        {"convert", [](){return std::unique_ptr<Command>(new Convert());}},
        {"density", [](){return std::unique_ptr<Command>(new Density());}},
        {"elastic", [](){return std::unique_ptr<Command>(new Elastic());}},
        {"formats", [](){return std::unique_ptr<Command>(new Formats());}},
        {"hbonds", [](){return std::unique_ptr<Command>(new HBonds());}},
        {"info", [](){return std::unique_ptr<Command>(new Info());}},
        {"merge", [](){return std::unique_ptr<Command>(new Merge());}},
        {"msd", [](){return std::unique_ptr<Command>(new MSD());}},
        {"rdf", [](){return std::unique_ptr<Command>(new Rdf());}},
        {"rotcf", [](){return std::unique_ptr<Command>(new Rotcf());}},
    };
    return commands;
}
