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
    /// `nx` and `ny`, and which can hold data in the `min_x - max_x` range (resp `min_y - max_y`).
    Histogram(size_t nx, double min_x, double max_x, size_t ny, double min_y, double max_y): 
        super(nx*ny),
        nx_(nx),
        min_x_(min_x),
        dx_((max_x - min_x) / nx),
        ny_(ny),
        min_y_(min_y),
        dy_((max_y - min_y) / ny) {
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
        return (*this)[j + i * ny_];
    }

    /// Get the size of the bins in the first dimension
    double dx() const {return dx_;}
    /// Get the size of the bins in the second dimension
    double dy() const {return dy_;}

    /// Get the minimal value of this histogram in the first dimension
    double min_x() const {return min_x_;}
    /// Get the minimal value of this histogram in the second dimension
    double min_y() const {return min_y_;}

    /// Get the number of bins in the first dimension
    double nx() const {return nx_;}
    /// Get the number of bins in the second dimension
    double ny() const {return ny_;}

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
        (*this)[bin_y + bin_x * ny_] += 1;
    }

    /// Get the x value corresponding to the ith element of the histogram
    T x(size_t i, bool radial = false) const {
        auto coord_x = min_x_ + i * dx_;
        if (radial and coord_x == 0) {
            coord_x = dx_ / 1000; // use a small r compared to dr to avoid Nan in the output
        }
        return coord_x;
    } 

    /// Get the y value corresponding to the ith element of the histogram
    T y(size_t i, bool radial = false) const {
        auto coord_y = min_y_ + i * dy_;
        if (radial and coord_y == 0) {
            coord_y = dy_ / 1000; // use a small r compared to dr to avoid Nan in the output
        }
        return coord_y;
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
    /// Number of bins in the first dimension
    size_t nx_ = 0;
    /// Width of a bin in the first dimension
    double dx_ = 0;
    /// Starting value for the histogram in the first dimension
    double min_x_ = 0;
    /// Number of bins in the second dimension
    size_t ny_ = 0;
    /// Width of a bin in the second dimension
    double dy_ = 0;
    /// Starting value for the histogram in the second dimension
    double min_y_ = 0;
};

#endif
