/* cfiles, an analysis frontend for the Chemfiles library
 * Copyright (C) 2015 Guillaume Fraux
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/
*/
#pragma once
#ifndef CFILES_RDF_HPP
#define CFILES_RDF_HPP

#include <vector>
#include <fstream>

#include "chemfiles.hpp"

#include "Histogram.hpp"
#include "Command.hpp"

struct rdf_options {
    //! Input trajectory
    std::string infile;
    //! Output data file
    std::string outfile;
    //! Selection for the atoms in radial distribution
    std::string selection;
    //! Number of points in the histogram
    size_t npoints;
    //! Maximum distance for the histogram
    double rmax;
    //! First step to use in g(r)
    size_t start;
    //! Last step to use in g(r)
    size_t end;
    //! Use a step every `stride` steps
    size_t stride;
    //! Unit cell to use
    std::vector<double> cell;
    //! Topology file to use
    std::string topology;
};

class Rdf : public Command {
public:
    Rdf(): selection_("all") {}
    virtual int run(int argc, const char* argv[]) override;
    virtual std::string description() const override;
    virtual std::string help() const override;
private:
    //! Add the data from a frame to the histogram
    void accumulate(const chemfiles::Frame& frame);
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
    //! Number of pairs found inside the rdf
    size_t nsteps_ = 0;

    //! Selection for the atoms in the pair
    chemfiles::Selection selection_;
};

#endif
