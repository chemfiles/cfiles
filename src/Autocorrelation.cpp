// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <numeric>
#include <cmath>
#include <cassert>

#include "Autocorrelation.hpp"

Autocorrelation::Autocorrelation(size_t size):
    size_(size),
#ifdef CFILES_USE_FFTW3
    fft_size_(2 * size_),
#else
    fft_size_(std::max(2 * size_, static_cast<size_t>(kiss_fftr_next_fast_size_real(size_)))),
#endif
    n_timeseries_(0),
    result_(size_, 0),
    spectrum_(fft_size_ / 2 + 1),
    direct_(fft_size_, false),
    reverse_(fft_size_, true)
{}

void Autocorrelation::add_timeserie(std::vector<float> timeserie) {
    // The algorithm used here compute autocorrelation using FFT.
    // It is described in https://doi.org/10.1016/0010-4655(95)00048-K
    assert(size_ == timeserie.size());
    n_timeseries_ += 1;

    // Pad the timeserie vector with 0 up to at least 2 * size_
    timeserie.insert(timeserie.end(), fft_size_ - size_, 0);

#ifdef CFILES_USE_FFTW3
    fftwf_execute_dft_r2c(direct_, timeserie.data(), spectrum_.data());
    for (auto& value: spectrum_) {
        value[0] = value[0] * value[0] + value[1] * value[1];
        value[1] = 0;
    }
    fftwf_execute_dft_c2r(reverse_, spectrum_.data(), timeserie.data());
#else
    kiss_fftr(direct_, timeserie.data(), spectrum_.data());
    for (auto& value: spectrum_) {
        // Replace values by their norm
        value.r = value.r * value.r + value.i * value.i;
        value.i = 0;
    }
    kiss_fftri(reverse_, spectrum_.data(), timeserie.data());
#endif

    for (size_t i=0; i<size_; i++) {
        result_[i] += timeserie[i];
    }
}
