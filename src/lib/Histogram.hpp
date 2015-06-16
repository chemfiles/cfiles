/*
 * chrp, an analysis frontend for the Chemharp library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#pragma once
#ifndef CHRP_FRONTEND_HISTOGRAM_HPP
#define CHRP_FRONTEND_HISTOGRAM_HPP

#include <cstdlib>
#include <cassert>
#include <vector>
#include <numeric>
#include <functional>
#include <type_traits>

#define ENABLE_FOR_ARITHMETIC_TYPES \
template <class U = T, typename = typename std::enable_if<std::is_arithmetic<U>::value>::type>

/** @class Histogram Histogram.hpp Histogram.cpp
 *
 * Histogram class, for averaging over a set of data. The data can be multi-dimmensional,
 * if T is defined to be std::array<U, N>, or any other array type. If the data is not
 * an arithmetic type (integer or floating point), then the position of the data to insert
 * must be specified by the caller. Else, the position of the data is automatically
 * determined by the bin_size and the
 */
template <class T>
class Histogram {
public:
    //! Default constructor
    Histogram() : Histogram(0) {}
    //! Constructor with a specific number of bins \c nbins
    Histogram(size_t nbins) : data_(nbins) {}
    //! Constructor with a specific number of bins \c nbins, and which can hold data in
    //! the \c min - \c max range.
    //! This can only be used if T is arithmetic (integer or floating point type).
    ENABLE_FOR_ARITHMETIC_TYPES
    Histogram(size_t nbins, double min, double max) : data_(nbins), dr_((max-min)/nbins) {}
    //! TODO
    //! This can only be used if T is arithmetic (integer or floating point type).
    ENABLE_FOR_ARITHMETIC_TYPES
    Histogram(size_t nbins, double dr) : data_(nbins), dr_(dr) {}

    ~Histogram() = default;
    Histogram(const Histogram&) = default;
    Histogram(Histogram&&) = default;
    Histogram& operator=(const Histogram&) = default;
    Histogram& operator=(Histogram&&) = default;

    //! Get the size of the bins
    ENABLE_FOR_ARITHMETIC_TYPES
    double bin_size() const {return dr_;}
    //! Set the size of the bins. This should be called before any attempt to add data
    //! inside the Histogram.
    ENABLE_FOR_ARITHMETIC_TYPES
    void bin_size(double dr) {dr_ = dr;}

    //! Insert some \c data in the histogram, at the position \c bin.
    void insert(T new_data, size_t bin) {
        assert(bin < data_.size());
        data_[bin] += new_data;
    }
    //! Insert some \c data in the histogram, and guess the position using the \c bin_size
    //! of the Histogram.
    ENABLE_FOR_ARITHMETIC_TYPES
    void insert(T new_data) {
        size_t bin = static_cast<size_t>(new_data / dr_);
        assert(bin < data_.size());
        data_[bin] += new_data;
    }

    //! Normalize the data so that the mean of the data is 1
    ENABLE_FOR_ARITHMETIC_TYPES
    void normalize() {
        double mean = std::accumulate(begin(data_), end(data_), 0.0) / data_.size();
        for (auto& val : data_) {
            val /= mean;
        }
    }
    //! Normalize the data with a \c function callback, which will be called for each value.
    //! The function should take two arguments being the current bin and the data, and
    //! return the new data.
    void normalize(std::function<T(size_t, T)> function) {
        for (size_t i=0; i< data_.size(); i++){
            data_[i] = function(i, data_[i]);
        }
    }

    //! Access the underlying data.
    const std::vector<T>& data() const {return data_;}
    //! Access the size of the underlying data.
    size_t size() const {return data_.size();}
    //! Remove all data.
    void clear() {data_.clear();}

    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    iterator begin() { return std::begin(data_);}
    const_iterator begin() const { return std::begin(data_);}
    iterator end() { return std::end(data_);}
    const_iterator end() const { return std::end(data_);}
private:
    //! Holding the data
    std::vector<T> data_;
    //! TODO
    double dr_ = 0;
    //! Number of elements added into this Histogram
    size_t nadded_ = 0;
};

#undef ENABLE_FOR_ARITHMETIC_TYPES

#endif
