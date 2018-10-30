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
        /// Selection/3D vector description for the optional center point
        std::string center;
        /// Number of points in the histogram
        size_t npoints = 0;
        /// Maximum distance for the histogram
        double rmax = 0;
    };

    Rdf(): selection_("all") {}
    std::string description() const override;

    Averager setup(int argc, const char* argv[]) override;
    void accumulate(const chemfiles::Frame& frame, Histogram& histogram) override;
    void finish(const Histogram& histogram) override;

private:
    /// Check if the maximal distance is larger than the biggest inscribed
    /// sphere in the frame unit cell
    void check_rmax(const chemfiles::Frame& frame) const;

    /// Options for this instance of RDF
    Options options_;
    /// Selection for the atoms in the pair
    chemfiles::Selection selection_;
    /// Selection for the center point
    chemfiles::optional<chemfiles::Selection> center_sel_ = chemfiles::nullopt;
    /// Fixed center point
    chemfiles::optional<chemfiles::Vector3D> center_ = chemfiles::nullopt;
    /// Also compute and average coordination numbers, for both i->j pairs and
    /// j->i pairs
    Averager coord_ij_;
    Averager coord_ji_;
};

#endif
