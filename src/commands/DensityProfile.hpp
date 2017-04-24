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
        /// Axis along which the profile should be computed
        std::string axis_str;
        size_t axis_vec[3];
        /// Number of points in the profile
        size_t npoints;
        /// Maximum in the profile
        double max = 0;
    };

    DensityProfile(): selection_("atoms: all"), axis_(0,0,1) {}
    std::string description() const override;

    Averager<double> setup(int argc, const char* argv[]) override;
    void accumulate(const chemfiles::Frame& frame, Histogram<double>& histogram);
    void finish(const Histogram<double>& histogram);

private:
    Options options_;
    chemfiles::Selection selection_;
    Axis axis_;
};

#endif
