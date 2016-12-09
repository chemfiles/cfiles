// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#ifndef CFILES_ERRORS_HPP
#define CFILES_ERRORS_HPP

#include <stdexcept>

struct CFilesError : public std::runtime_error {
    CFilesError(const std::string& message): std::runtime_error(message) {}
};

#endif
