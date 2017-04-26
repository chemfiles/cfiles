// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_AXIS_HPP
#define CFILES_AXIS_HPP

#include <chemfiles.hpp>
#include "Errors.hpp"
#include <cmath>
#include <iostream>

using namespace chemfiles;

class Axis {
public:
    /// Constructor with an axis name (X,Y,Z, XY, XZ, YZ...)
    explicit Axis(std::string name): Axis(0,0,0) {
        if (name=="X" or name=="x") {
            vector_ = {1,0,0};
        } else if (name=="Y" or name=="y") {
            vector_ = {0,1,0};
        } else if (name=="Z" or name=="z") {
            vector_ = {0,0,1};
        } else {
            throw CFilesError("Axis non implemented, enter vector coordinates instead");
        }
    }
    /// Constructor with a 3D vector
    Axis(double a, double b, double c): vector_(vector3d(a,b,c)) {
        /// normalize the axis
        vector_ = vector_ / norm(vector_);
    }

    /// print axis
    Vector3D& get_coordinates() { return vector_; }
    bool is_null() { return vector_ == vector3d(0,0,0); }

    /// projection on axis (may be negative)
    double projection(const Vector3D & positions) {
        auto dot_product = dot(vector_,positions);
        return dot_product;
    }

    /// radial distance to axis
    double radial(const Vector3D & positions) {
        auto dot_product = dot(vector_,positions);
        auto norme2 = norm2(positions);
        auto distance2 = norme2 - dot_product * dot_product;
        return sqrt(distance2);
    }

private:
    /// axis coordinates
    Vector3D vector_;
};

#endif
