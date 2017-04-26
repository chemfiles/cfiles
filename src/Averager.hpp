// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_AVERAGER_HPP
#define CFILES_AVERAGER_HPP

#include "Histogram.hpp"

/// Average class, averaging an historgram over multiple steps
template <class T>
class Averager: public Histogram<T> {
public:
    /// Default constructor
    Averager(): Histogram<T>(), averaged_() {}
    /// Constructor with a specific number of bins `nbins`, and which can hold
    /// data in the `min - max` range.
    Averager(size_t nbins, double min, double max): Histogram<T>(nbins, min, max), averaged_(nbins) {}

    Averager(const Averager&) = default;
    Averager(Averager&&) = default;
    Averager& operator=(const Averager&) = default;
    Averager& operator=(Averager&&) = default;

    /// Store the current data for averaging, and clean the current data
    /// (set it to `T()`)
    void step() {
        for (size_t i=0; i<this->size(); i++) {
            averaged_[i] += (*this)[i];
            (*this)[i] = T();
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
    std::vector<T> averaged_;
    /// Number of time `step` was called
    size_t nsteps_ = 0;
};

#endif
