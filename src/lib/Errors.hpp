/*
 * chrp, an analysis frontend for the Chemharp library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#pragma once
#ifndef CHRP_FRONTEND_ERRORS_HPP
#define CHRP_FRONTEND_ERRORS_HPP

#include <stdexcept>

class chrp_exception : public std::runtime_error {
    using runtime_error::runtime_error;
};

#endif
