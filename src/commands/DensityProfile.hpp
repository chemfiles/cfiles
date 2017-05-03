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
        /// Type of profile: 0 means none, 1 linear profile, 2 radial profile
        size_t type_profile[2] = {0, 0};
	/// Coordinate of origin
        Vector3D origin[2];
        /// Number of points in the profile
        size_t npoints[2];
        /// Maximum in the profile
        double max[2] = {0, 0};
        /// Minimum in the profile
        double min[2] = {0, 0};
    };

    DensityProfile(): selection_("atoms: all"), axis_x_(0,0,1), axis_y_(0,0,1) {}
    std::string description() const override;

    Averager<double> setup(int argc, const char* argv[]) override;
    void accumulate(const chemfiles::Frame& frame, Histogram<double>& histogram);
    void finish(const Histogram<double>& histogram);

private:
    Options options_;
    chemfiles::Selection selection_;
    size_t n_axis_;
    Axis axis_x_;
    Axis axis_y_;
};

#endif
