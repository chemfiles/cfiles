// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#ifndef CFILES_HISTOGRAM_HPP
#define CFILES_HISTOGRAM_HPP

#include <vector>
#include <cmath>
#include <numeric>
#include <functional>

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
    Histogram(): Histogram(0, 0) {}
    /// Constructor with a specific number of bins `nbins`, and which can hold
    /// data in the `min - max` range.
    Histogram(size_t nbins, double min, double max): Histogram(nbins, (max - min) / nbins) {}
    /// Constructor with a specific number of `bins` and a specific bin size `dr`
    Histogram(size_t nbins, double dr): super(nbins), dr_(dr) {
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

    /// Insert some `data` in the histogram, and guess the position using the
    /// `bin_size` of the Histogram.
    void insert(T new_data) {
        auto bin = std::floor(new_data / dr_);
        assert(bin < size());
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
};

#endif
