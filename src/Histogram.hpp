// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_HISTOGRAM_HPP
#define CFILES_HISTOGRAM_HPP

#include <vector>
#include <cmath>
#include <numeric>
#include <functional>

#include <fmt/format.h>
#include "warnings.hpp"

/// Histogram class
class Histogram {
public:
    using iterator = std::vector<double>::const_iterator;

    /// Information for each dimension of a Histogram
    struct Dimension {
        Dimension(size_t n, double min, double max): nbins(n), start(min), width((max - min) / n) {}

        /// Number of bins
        size_t nbins;
        /// Starting value for the histogram
        double start;
        /// Width of a bin
        double width;

        double stop() const {
            return start + nbins * width;
        }

        double coord(size_t i) const {
            return start + (i + 0.5) * width;
        }
    };

    /// Default constructor
    Histogram(): Histogram(0, 0, 0, 0, 0, 0) {}
    /// Constructor for a flat 2d histogram with a specific number of bins in
    /// each direction `n1` and `n2`, and which can hold data in the `min_1 -
    /// mafirst_2` range (resp `min_1 - mafirst_2`).
    Histogram(size_t n1, double min1, double max1, size_t n2, double min2, double max2)
        : data_(n1 * n2), first_(n1, min1, max1), second_(n2, min2, max2) {}

    /// Constructor for a 1d histogram with a specific number of bins `n_bins`,
    /// and which can hold data in the `min - max` range.
    Histogram(size_t n_bins, double min, double max): Histogram(n_bins, min, max, 1, 0, 1) {}

    Histogram(const Histogram&) = default;
    Histogram(Histogram&&) = default;
    Histogram& operator=(const Histogram&) = default;
    Histogram& operator=(Histogram&&) = default;

    size_t size() const {
        return data_.size();
    }

    iterator begin() const {
        return data_.begin();
    }

    iterator end() const {
        return data_.end();
    }

    double operator[](size_t i) const {
        return data_[i];
    }

    double& operator[](size_t i) {
        return data_[i];
    }

    /// Using call pperator for 2D indexing 2D histogram
    double operator()(size_t i, size_t j) const {
        return data_[j + i * second_.nbins];
    }

    /// Get the first dimension
    const Dimension& first() const {return first_;}

    /// Get the second dimension
    const Dimension& second() const {return second_;}

    /// Insert some `x,y` in the histogram
    void insert(double x, double y = 0) {
        auto bin1 = std::floor((x - first_.start) / first_.width);
        auto bin2 = std::floor((y - second_.start) / second_.width);
        if (bin1 >= first_.nbins or bin1 < 0) {
            warn_once(fmt::format(
                "point {} is out of histogram boundaries ({}:{})",
                x, first_.start, first_.stop()
            ));
            return;
        }
        if (bin2 >= second_.nbins or bin2 < 0) {
            warn_once(fmt::format(
                "point {} is out of histogram boundaries ({}:{})",
                y, second_.start, second_.stop()
            ));
            return;
        }
        data_[bin2 + bin1 * second_.nbins] += 1;
    }

    /// Normalize the data with a `function` callback, which will be called for
    /// each value. The function should take two arguments being the current
    /// bin index and the data, and return the new data.
    void normalize(std::function<double(size_t, double)> function) {
        for (size_t i = 0; i < this->size(); i++){
            data_[i] = function(i, data_[i]);
        }
    }
private:
    /// Histogram data
    std::vector<double> data_;
    /// First dimension
    Dimension first_;
    /// Second dimension
    Dimension second_;
};

#endif
