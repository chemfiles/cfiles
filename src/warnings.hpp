// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_WARNINGS_HPP
#define CFILES_WARNINGS_HPP

#include <string>

/// Print a warning to the standard error stream
void warn(std::string message);

/// Print a warning once to the standard error stream
void warn_once(std::string message);

#endif
