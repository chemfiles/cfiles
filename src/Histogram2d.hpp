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
    Histogram(): Histogram(0, 0, 0, 0, 0, 0) {}
    /// Constructor with a specific number of bins `nbins`, and which can hold
    /// data in the `min - max` range.
    Histogram(size_t nbins1, double min1, double max1, size_t nbins2, double min2, double max2): super(nbins1,nbins2), min1_(min1), min2_(min2), dr1_((max1 - min1) / nbins1), dr2_((max2 - min2) / nbins2) {
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
    double bin_size1() const {return dr1_;}
    double bin_size2() const {return dr2_;}

    /// Get the min
    double min1() const {return min1_;}
    double min2() const {return min2_;}

    /// Insert some `data` in the histogram, and guess the position using the
    /// `bin_size` of the Histogram.
    void insert(T new_data) {
        auto bin = std::floor((new_data - min_) / dr_);
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
    double dr1_ = 0;
    double dr2_ = 0;
    double min1_ = 0;
    double min2_ = 0;
};

#endif
