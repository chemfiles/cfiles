// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_DEFAULT_OPTIONS_HPP
#define CFILES_DEFAULT_OPTIONS_HPP

#include <Options.hpp>

Option OPTION_OUTPUT(const std::string& extension);

extern Option OPTION_FORMAT;
extern Option OPTION_TOPOLOGY;
extern Option OPTION_TOPOLOGY_FORMAT;
extern Option OPTION_GUESS_BONDS;
extern Option OPTION_CELL;
extern Option OPTION_STEPS;

#endif
