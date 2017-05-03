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
    /// Constructor for a flat 2d histogram with a specific number of bins in each direction 
    /// `n_x` and `n_y`, and which can hold data in the `min_x - max_x` range (resp `min_y - max_y`).
    Histogram(size_t n_x, double min_x, double max_x, size_t n_y, double min_y, double max_y): super(n_x*n_y) {
        nx_ = n_x;
        min_x_ = min_x;
        dx_ = (max_x - min_x) / n_x;
        ny_ = n_y;
        min_y_ = min_y;
        dy_ = (max_y - min_y) / n_y;
        static_assert(
            std::is_arithmetic<T>::value,
            "Histogram<T> is only defined for arithmetics types T"
        );
    }

    /// Constructor for a 1d histogram with a specific number of bins `n_bins`, and which can hold 
    /// data in the `min - max` range.
    Histogram(size_t n_bins, double min, double max): Histogram(n_bins, min, max, 1, 0, 1) {}

    Histogram(const Histogram&) = default;
    Histogram(Histogram&&) = default;
    Histogram& operator=(const Histogram&) = default;
    Histogram& operator=(Histogram&&) = default;

    /// Get the size of the bins
    double dx() const {return dx_;}
    double dy() const {return dy_;}

    /// Get the minimal value of this histogram
    double min_x() const {return min_x_;}
    double min_y() const {return min_y_;}

    /// Insert some `x,y` in the histogram
    void insert(T x, T y = 0) {
        auto bin_x = std::floor((x - min_x_) / dx_);
        auto bin_y = std::floor((y - min_y_) / dy_);
        if (bin_x >= nx_ or bin_x < 0) {
            std::string s = std::to_string(x);
            throw OutOfBoundsError("Element " + s + " out of boundaries");
        }
        if (bin_y >= ny_ or bin_y < 0) {
            std::string s = std::to_string(y);
            throw OutOfBoundsError("Element " + s + " out of boundaries");
        }
        (*this)[bin_x + bin_y * nx_] += 1;
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
    size_t nx_ = 0;
    double dx_ = 0;
    double min_x_ = 0;
    size_t ny_ = 0;
    double dy_ = 0;
    double min_y_ = 0;
};

#endif
