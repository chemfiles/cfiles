/* cfiles, an analysis frontend for the Chemfiles library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
#pragma once
#ifndef CFILES_HELP_HPP
#define CFILES_HELP_HPP

#include "Command.hpp"

class Help : public Command {
public:
    Help() = default;
    virtual int run(int argc, char** argv) override;
    virtual std::string description() const override;
private:
    //! List all available commands with the associated description
    void list_commands() const;
    //! Get help about one command
    void about(const std::string& command) const;
};

#endif
