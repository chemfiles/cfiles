// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_DENSITY_PROFILE_HPP
#define CFILES_DENSITY_PROFILE_HPP

#include <chemfiles.hpp>

#include "AveCommand.hpp"
#include "Axis.hpp"
#include "utils.hpp"

class DensityProfile final: public AveCommand {
public:
    struct Options {
    /// Output
    std::string outfile;
    /// Selection for the donor-acceptor
    std::string selection;
    /// Coordinate of origin
        Vector3D origin;
        /// Number of points in the profile
        size_t npoints[2];
        /// Maximum in the profile
        double max[2] = {0, 0};
        /// Minimum in the profile
        double min[2] = {0, 0};
    };

    DensityProfile(): selection_("atoms: all"), axis_() {}
    std::string description() const override;

    Averager setup(int argc, const char* argv[]) override;
    void accumulate(const chemfiles::Frame& frame, Histogram& histogram) override;
    void finish(const Histogram& histogram) override;

    size_t dimensionality() { return axis_.size();}

private:
    Options options_;
    chemfiles::Selection selection_;
    std::vector<Axis> axis_;
};

#endif
