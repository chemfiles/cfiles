/* cfiles, an analysis frontend for the Chemfiles library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#pragma once
#ifndef CFILES_CONVERT_HPP
#define CFILES_CONVERT_HPP

#include "Command.hpp"

class Convert : public Command {
public:
    Convert() {}
    virtual int run(int argc, char** argv) override;
    virtual std::string description() const override;
    virtual std::string help() const override;
private:
};

#endif
