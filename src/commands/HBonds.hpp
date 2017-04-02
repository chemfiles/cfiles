// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#ifndef CFILES_HBONDS_HPP
#define CFILES_HBONDS_HPP

#include <chemfiles.hpp>

#include "Command.hpp"
#include "utils.hpp"

class HBonds final: public Command {
public:
    struct Options {
        /// Input trajectory
        std::string trajectory;
	/// Output
	std::string outfile;
        /// Specific format to use with the trajectory
        std::string format = "";
        /// Specific steps to use from the trajectory
        steps_range steps;
	/// Selection for the donor-acceptor
	std::string selection;
        /// Do we have a custom cell to use?
        bool custom_cell = false;
        /// Unit cell to use
        chemfiles::UnitCell cell;
        /// Topology file to use
        std::string topology = "";
        /// Format to use for the topology file
        std::string topology_format = "";
        /// Should we try to guess the topology?
        bool guess_bonds = false;
	/// Parameters for donor-acceptor max distance (in angstroms)
	/// and for donor-acceptor-hydrogen max angle (in degrees)
	double parameters[2] = {3.0, 30.0};
    };

    HBonds(): selectionAcceptor_("bonds: type(#2) == H"), selectionDonor_("atoms: all") {}
    int run(int argc, const char* argv[]) override;
    std::string description() const override;

private:
    Options options_;
    chemfiles::Selection selectionAcceptor_;
    chemfiles::Selection selectionDonor_;
};

#endif
