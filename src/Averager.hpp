// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_AVERAGER_HPP
#define CFILES_AVERAGER_HPP

#include "Histogram.hpp"

/// Average class, averaging an historgram over multiple steps
class Averager: public Histogram {
public:
    /// Default constructor
    Averager(): Histogram(), averaged_() {}
    /// Constructor for a flat 2d histogram with a specific number of bins in each direction
    /// `n1` and `n2`, and which can hold data in the `min1 - max1` range (resp `min2 - max2`).
    Averager(size_t n1, double min1, double max1, size_t n2, double min2, double max2):
        Histogram(n1, min1, max1, n2, min2, max2), averaged_(n1 * n2) {}
    /// Constructor with a specific number of bins `nbins`, and which can hold
    /// data in the `min - max` range.
    Averager(size_t nbins, double min, double max): Histogram(nbins, min, max), averaged_(nbins) {}

    Averager(const Averager&) = default;
    Averager(Averager&&) = default;
    Averager& operator=(const Averager&) = default;
    Averager& operator=(Averager&&) = default;

    /// Store the current data for averaging, and clean the current data
    /// (set it to `T()`)
    void step() {
        for (size_t i=0; i<this->size(); i++) {
            averaged_[i] += (*this)[i];
            (*this)[i] = 0;
        }
        nsteps_++;
    }

    void average() {
        for (size_t i=0; i<this->size(); i++) {
            (*this)[i] = averaged_[i] / nsteps_;
        }
    }

private:
    /// Accumulating the averaged values
    std::vector<double> averaged_;
    /// Number of time `step` was called
    size_t nsteps_ = 0;
};

#endif
