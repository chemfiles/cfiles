// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_AUTOCORRELATION_HPP
#define CFILES_AUTOCORRELATION_HPP

#include <kiss_fftr.h>
#include <vector>

#include "Errors.hpp"

/// A RAII capsule for kiss_fftr_cfg
class fft_config {
public:
    fft_config():fft_config(0, false) {}

    fft_config(size_t size, bool reverse) {
        cfg_ = kiss_fftr_alloc(size, reverse, nullptr, nullptr);
        if (cfg_ == nullptr) {
            throw CFilesError("Could not alocate memory for FFT");
        }
    }

    ~fft_config() {
        kiss_fft_free(cfg_);
    }

    fft_config(fft_config&& other): cfg_(other.cfg_) {
        other.cfg_ = nullptr;
    }

    fft_config& operator=(fft_config&& other) {
        kiss_fft_free(this->cfg_);
        this->cfg_ = other.cfg_;
        other.cfg_ = nullptr;
        return *this;
    }

    operator kiss_fftr_cfg() const {
        return cfg_;
    }
private:
    kiss_fftr_cfg cfg_ = nullptr;
};


class Autocorrelation {
public:
    Autocorrelation(size_t size);

    /// Compute autocorrelation for the given time serie, and store it for
    /// future averaging
    void add_timeserie(std::vector<float> timeserie);

    /// Normalize the averaged autocorrelations
    void normalize() {
        for (size_t i=0; i<size_; i++) {
            result_[i] /= n_timeseries_;
        }
    }

    /// Get the averaged autocorrelations
    const std::vector<float>& average() const {
        return result_;
    }

private:
    /// Compute autocorrelation using FFT algorithm, and store the result in
    /// `timeserie`
    void compute_fft(std::vector<float>& timeserie);

    /// Compute autocorrelation using direct algorithm, and store the result in
    /// `timeserie`
    void compute_direct(std::vector<float>& timeserie);

    /// Number of elements in the time series
    size_t size_;
    /// Number of timeseries used
    size_t n_timeseries_ = 0;
    /// Accumulated autocorrelations
    std::vector<float> result_;

    /// Buffer for FFT data
    std::vector<kiss_fft_cpx> buffer_;
    /// FFT configuration
    fft_config direct_;
    /// Reverse FFT configuration
    fft_config reverse_;
};

#endif
