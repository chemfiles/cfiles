// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_AUTOCORRELATION_HPP
#define CFILES_AUTOCORRELATION_HPP

#include <vector>
#include "Errors.hpp"

#ifdef CFILES_USE_FFTW3
#include <fftw3.h>
#else
#include <kiss_fftr.h>
#endif

#ifdef CFILES_USE_FFTW3
/// A RAII capsule for fftwf_plan
class FFTWPlan {
public:
    FFTWPlan(size_t size, bool reverse) {
        if (reverse) {
            // Use the FFTW_UNALIGNED flag, as this will be used for multiple
            // array, over which this code does not have control w.r.t.
            // alignment. Allocating memory with fftw_malloc would need a big
            // refactoring.
            plan_ = fftwf_plan_dft_c2r_1d(size, nullptr, nullptr, FFTW_ESTIMATE | FFTW_UNALIGNED);
        } else {
            plan_ = fftwf_plan_dft_r2c_1d(size, nullptr, nullptr, FFTW_ESTIMATE | FFTW_UNALIGNED);
        }
        if (plan_ == nullptr) {
            throw CFilesError("Could not allocate memory for FFT");
        }
    }

    ~FFTWPlan() {
        fftwf_destroy_plan(plan_);
    }

    FFTWPlan(FFTWPlan&& other): plan_(other.plan_) {
        other.plan_ = nullptr;
    }

    FFTWPlan& operator=(FFTWPlan&& other) {
        fftwf_destroy_plan(this->plan_);
        this->plan_ = other.plan_;
        other.plan_ = nullptr;
        return *this;
    }

    operator fftwf_plan() const {
        return plan_;
    }
private:
    fftwf_plan plan_ = nullptr;
};

using fft_complex = fftwf_complex;
using fft_plan = FFTWPlan;

#else
/// A RAII capsule for kiss_fftr_cfg
class KissFTTConfig {
public:
    KissFTTConfig(size_t size, bool reverse) {
        cfg_ = kiss_fftr_alloc(size, reverse, nullptr, nullptr);
        if (cfg_ == nullptr) {
            throw CFilesError("Could not allocate memory for FFT");
        }
    }

    ~KissFTTConfig() {
        kiss_fft_free(cfg_);
    }

    KissFTTConfig(KissFTTConfig&& other): cfg_(other.cfg_) {
        other.cfg_ = nullptr;
    }

    KissFTTConfig& operator=(KissFTTConfig&& other) {
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

using fft_complex = kiss_fft_cpx;
using fft_plan = KissFTTConfig;
#endif

class Autocorrelation {
public:
    Autocorrelation(size_t size);
    ~Autocorrelation();

    Autocorrelation(const Autocorrelation&) = delete;
    Autocorrelation& operator=(const Autocorrelation&) = delete;

    Autocorrelation(Autocorrelation&&);
    Autocorrelation& operator=(Autocorrelation&&);

    /// Compute autocorrelation for the given time serie, and store it for
    /// future averaging
    void add_timeserie(std::vector<float> timeserie);

    /// Normalize the averaged autocorrelations
    void normalize() {
        for (size_t i=0; i<size_; i++) {
            // fft_size_ is the gain from doing FFT -> iFFT with both FFTW3 and
            // KissFFT
            result_[i] /=  fft_size_ * n_timeseries_ * (size_ - i);
        }
    }

    /// Get the averaged autocorrelations
    const std::vector<float>& get_result() const {
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
    /// Number of points for the FFT
    size_t fft_size_;
    /// Number of timeseries used
    size_t n_timeseries_;
    /// Accumulated autocorrelations
    std::vector<float> result_;

    /// Buffer for FFT data
    fft_complex* spectrum_;
    /// FFT configuration
    fft_plan direct_;
    /// Reverse FFT configuration
    fft_plan reverse_;
};

#endif
