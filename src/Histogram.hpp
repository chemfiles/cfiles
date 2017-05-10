// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_HISTOGRAM_HPP
#define CFILES_HISTOGRAM_HPP

#include <vector>
#include <cmath>
#include <numeric>
#include <functional>
#include "Errors.hpp"

/// Information for each dimension of a Histogram
struct HistDim {
    /// Number of bins
    size_t nbins;
    /// Starting value for the histogram
    double min;
    /// Width of a bin
    double dr;
};

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
    /// `n1` and `n2`, and which can hold data in the `min_1 - max_2` range (resp `min_1 - max_2`).
    Histogram(size_t n1, double min1, double max1, size_t n2, double min2, double max2): 
        super(n1*n2),
        first_dimension_{n1, min1, (max1 - min1) / n1},
        second_dimension_{n2, min2, (max2 - min2) / n2} {
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
    /// Operator to get the (i,j) element of a 2D histogram
    const T& operator()(size_t i, size_t j) const {
        return (*this)[j + i * second_dimension_.nbins];
    }

    /// Get the first dimension
    const HistDim& first() const {return first_dimension_;}

    /// Get the second dimension
    const HistDim& second() const {return second_dimension_;}

    /// Insert some `x,y` in the histogram
    void insert(T x, T y = 0) {
        auto bin1 = std::floor((x - first_dimension_.min) / first_dimension_.dr);
        auto bin2 = std::floor((y - second_dimension_.min) / second_dimension_.dr);
        if (bin1 >= first_dimension_.nbins or bin1 < 0) {
            std::string s = std::to_string(x);
            throw OutOfBoundsError("Element " + s + " out of boundaries");
        }
        if (bin2 >= second_dimension_.nbins or bin2 < 0) {
            std::string s = std::to_string(y);
            throw OutOfBoundsError("Element " + s + " out of boundaries");
        }
        (*this)[bin2 + bin1 * second_dimension_.nbins] += 1;
    }

    /// Get the x value corresponding to the ith element of the histogram
    T first_index(size_t i, bool radial = false) const {
        auto coord = first_dimension_.min + i * first_dimension_.dr;
        if (radial and coord == 0) {
            coord = first_dimension_.min + (i + 0.5) * first_dimension_.dr; // shift the histogram to avoid Nan in the output
        }
        return coord;
    } 

    /// Get the y value corresponding to the ith element of the histogram
    T second_index(size_t i, bool radial = false) const {
        auto coord = second_dimension_.min + i * second_dimension_.dr;
        if (radial and coord == 0) {
            coord = second_dimension_.min + (i + 0.5) * second_dimension_.dr; // shift the histogram to avoid Nan in the output
        }
        return coord;
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
    /// First dimension
    HistDim first_dimension_;
    /// Second dimension
    HistDim second_dimension_;
};

#endif
