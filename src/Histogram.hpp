// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_HISTOGRAM_HPP
#define CFILES_HISTOGRAM_HPP

#include <vector>
#include <cmath>
#include <numeric>
#include <functional>
#include "Errors.hpp"

/// Histogram class
template <class T>
class Histogram: private std::vector<T> {
    using super = std::vector<T>;
public:
    using super::size;
    using super::operator[];
    using super::begin;
    using super::end;

    /// Default constructor
    Histogram(): Histogram(0, 0, 0) {}
    /// Constructor with a specific number of bins `nbins`, and which can hold
    /// data in the `min - max` range.
    Histogram(size_t nbins, double min, double max): super(nbins), min_(min), dr_((max - min) / nbins) {
        static_assert(
            std::is_arithmetic<T>::value,
            "Histogram<T> is only defined for arithmetics types T"
        );
    }

    Histogram(const Histogram&) = default;
    Histogram(Histogram&&) = default;
    Histogram& operator=(const Histogram&) = default;
    Histogram& operator=(Histogram&&) = default;

    /// Get the size of the bins
    double bin_size() const {return dr_;}

    /// Get the minimal value of this histogram
    double min() const {return min_;}

    /// Insert some `data` in the histogram
    void insert(T data) {
        auto bin = std::floor((data - min_) / dr_);
        if (bin >= size() or bin < 0) {
            throw OutOfBoundsError("Element out of boundaries");
        }
        (*this)[bin] += 1;
    }

    /// Normalize the data with a `function` callback, which will be called for
    /// each value. The function should take two arguments being the current
    /// bin index and the data, and return the new data.
    void normalize(std::function<T(size_t, T)> function) {
        for (size_t i=0; i< size(); i++){
            (*this)[i] = function(i, (*this)[i]);
        }
    }
private:
    /// Width of a bin in the Histogram, for arithmetic types only
    double dr_ = 0;
    double min_ = 0;
};

#endif
