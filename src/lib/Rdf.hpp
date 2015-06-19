/*
 * chrp, an analysis frontend for the Chemharp library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/

#pragma once
#ifndef CHRP_FRONTEND_RDF_HPP
#define CHRP_FRONTEND_RDF_HPP

#include <vector>

#include "Histogram.hpp"
#include "Command.hpp"

namespace harp {
    class Frame;
}

struct rdf_options {
    std::string infile;
    std::string outfile;
    size_t nbins;
    size_t start;
    size_t end;
    size_t stride;
    std::vector<double> cell;
    std::string topology;
};

class Rdf : public Command {
public:
    Rdf() : nsteps_(0) {}
    virtual int run(int argc, char** argv) override;
    virtual std::string description() override;
    virtual std::string help() override;
private:
    //! Add the data from a frame to the histogram
    void accumulate(harp::Frame& frame);
    //! Normalize the histogram data
    void finish();
    //! Write the histogram data to a file
    void write(const std::string& filename);

    //! Options for this instance of RDF
    rdf_options options_;
    //! Histogram, for inserting the data at each step
    Histogram<double> histogram_;
    //! Result for storing the pre-normalized results
    std::vector<double> result_;
    //! Number of steps we performed
    size_t nsteps_;
};

#endif
