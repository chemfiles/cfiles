// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_COMMANDS_HPP
#define CFILES_COMMANDS_HPP

#include <chemfiles.hpp>

#include "Command.hpp"
#include "Averager.hpp"

class Axis;

/// Base class for time-averaged computations
class AveCommand: public Command {
public:
    AveCommand(Options options);
    ~AveCommand() override = default;
    int run() override final;

    /// Setup the command and the histogram.
    virtual Averager setup() = 0;
    /// Add the data from a `frame` to the `histogram`
    virtual void accumulate(const chemfiles::Frame& frame, Histogram& histogram) = 0;
    /// Finish the run, and write any output
    virtual void finish(const Histogram& histogram) = 0;

protected:
    /// Trajectory path
    std::string trajectory_;

private:
    /// Averaging histogram for the data
    Averager histogram_;
};

/// Convert between file formats
class Convert final: public Command {
public:
    Convert();
    ~Convert() override = default;
    int run() override;
};

/// Compute elastic properties
class Elastic final: public Command {
public:
    Elastic();
    ~Elastic() override = default;
    int run() override;
};

/// Compute angles distributions
class AngleDistribution final: public AveCommand {
public:
    AngleDistribution();
    ~AngleDistribution() override = default;

    Averager setup() override;
    void accumulate(const chemfiles::Frame& frame, Histogram& histogram) override;
    void finish(const Histogram& histogram) override;

private:
    chemfiles::Selection selection_ = chemfiles::Selection("none");
};


class DensityProfile final: public AveCommand {
public:
    DensityProfile();
    ~DensityProfile() override = default;

    Averager setup() override;
    void accumulate(const chemfiles::Frame& frame, Histogram& histogram) override;
    void finish(const Histogram& histogram) override;

private:
    size_t dimensionality() const;

    chemfiles::Selection selection_ = chemfiles::Selection("all");
    chemfiles::Vector3D origin_ = chemfiles::Vector3D(0, 0, 0);
    std::vector<Axis> axis_;
};

#endif
