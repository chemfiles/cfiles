// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#ifndef CFILES_ANGLES_HPP
#define CFILES_ANGLES_HPP

#include "AveCommand.hpp"
#include "utils.hpp"

class AngleDistribution final: public AveCommand {
public:
    struct Options {
        /// Output data file
        std::string outfile;
        /// Selection for the atoms in radial distribution
        std::string selection;
        /// Number of points in the histogram
        size_t npoints;
    };

    AngleDistribution(): selection_("angles: all") {}
    std::string description() const override;

    Averager<double> setup(int argc, const char* argv[]) override;
    void accumulate(const chemfiles::Frame& frame, Histogram<double>& histogram) override;
    void finish(const Histogram<double>& histogram) override;

private:
    /// Options for this instance of RDF
    Options options_;
    /// Selection for the atoms in the pair
    chemfiles::Selection selection_;
};

#endif
