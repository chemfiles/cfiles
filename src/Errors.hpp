// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_ERRORS_HPP
#define CFILES_ERRORS_HPP

#include <stdexcept>
#include <fmt/format.h>

struct CFilesError : public std::runtime_error {
    CFilesError(const std::string& message): std::runtime_error(message) {}
};

struct OptionError : public CFilesError {
    OptionError(const std::string& message): CFilesError(message) {}
};

template <typename... Args>
inline CFilesError cfiles_error(const char *format, const Args & ... arguments) {
    return CFilesError(fmt::format(format, arguments...));
}

template <typename... Args>
inline OptionError option_error(const char *format, const Args & ... arguments) {
    return OptionError(fmt::format(format, arguments...));
}

#endif
