// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_ERRORS_HPP
#define CFILES_ERRORS_HPP

#include <stdexcept>

struct CFilesError : public std::runtime_error {
    CFilesError(const std::string& message): std::runtime_error(message) {}
};

#endif
