// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_AXIS_HPP
#define CFILES_AXIS_HPP

#include <chemfiles.hpp>
#include "Errors.hpp"
#include <cmath>

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
        if (vector_  != vector3d(0,0,0)) {
            vector_ = vector_ / norm(vector_);
        }
    }

    /// Get the vector coordinates
    Vector3D& get_coordinates() { return vector_; }
    /// Check if axis is null
    bool is_null() { return vector_ == vector3d(0,0,0); }

    /// projection on axis (may be negative)
    double projection(const Vector3D & positions) {
        return dot(vector_,positions);
    }

    /// radial distance to axis
    double radial(const Vector3D & positions) {
        auto dot_product = dot(vector_,positions);
        auto distance2 = norm2(positions) - dot_product * dot_product;
        return sqrt(distance2);
    }

private:
    /// axis coordinates
    Vector3D vector_;
};

#endif
