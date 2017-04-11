// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_RDF_HPP
#define CFILES_RDF_HPP

#include "AveCommand.hpp"

class Rdf final: public AveCommand {
public:
    struct Options {
        /// Output data file
        std::string outfile;
        /// Selection for the atoms in radial distribution
        std::string selection;
        /// Number of points in the histogram
        size_t npoints;
        /// Maximum distance for the histogram
        double rmax;
    };

    Rdf(): selection_("all") {}
    std::string description() const override;

    Averager<double> setup(int argc, const char* argv[]) override;
    void accumulate(const chemfiles::Frame& frame, Histogram<double>& histogram) override;
    void finish(const Histogram<double>& histogram) override;

private:
    /// Options for this instance of RDF
    Options options_;
    /// Selection for the atoms in the pair
    chemfiles::Selection selection_;
    /// Also compute and average coordination numbers
    Averager<double> coordination_;
};

#endif
