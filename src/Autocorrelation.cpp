// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <numeric>
#include <cmath>
#include <cassert>

#include "Autocorrelation.hpp"

Autocorrelation::Autocorrelation(size_t size):
    size_(size),
    result_(size_, 0),
    buffer_(2 * size_, {0}),
    direct_(2 * size_, false),
    reverse_(2 * size_, true) {}

void Autocorrelation::add_timeserie(std::vector<float> timeserie) {
    // The algorithm used here compute autocorrelation using FFT.
    // It is described in https://doi.org/10.1016/0010-4655(95)00048-K
    assert(size_ == timeserie.size());
    n_timeseries_ += 1;

    // Pad the timeserie vector with 0 to 2N
    timeserie.insert(timeserie.end(), size_, 0);

    kiss_fftr(direct_, timeserie.data(), buffer_.data());
    // Replace values by their norm
    for (auto& value: buffer_) {
        value.r = value.r * value.r + value.i * value.i;
        value.i = 0;
    }
    kiss_fftri(reverse_, buffer_.data(), timeserie.data());

    float norm = timeserie[0] / size_;
    for (size_t i=0; i<size_; i++) {
        result_[i] += timeserie[i] / (norm * (size_ - i));
    }
}
