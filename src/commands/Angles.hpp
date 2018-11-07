// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_ANGLES_HPP
#define CFILES_ANGLES_HPP

#include "AveCommand.hpp"
#include "utils.hpp"

class Angles final: public AveCommand {
public:
    struct Options {
        /// Output data file
        std::string outfile;
        /// Selection for the atoms in radial distribution
        std::string selection;
        /// Number of points in the histogram
        size_t npoints;
    };

    Angles(): selection_("angles: all") {}
    std::string description() const override;

    Averager setup(int argc, const char* argv[]) override;
    void accumulate(const chemfiles::Frame& frame, Histogram& histogram) override;
    void finish(const Histogram& histogram) override;

private:
    /// Options for this instance of RDF
    Options options_;
    /// Selection for the atoms in the pair
    chemfiles::Selection selection_;
};

#endif
