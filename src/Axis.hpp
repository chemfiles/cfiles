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
        if (name=="X") {
            vector_[0] = 1;
        } else if (name=="Y") {
            vector_[1] = 1;
        } else if (name=="Z") {
            vector_[2] = 1;
        } else {
            throw CFilesError("Axis non implemented, enter vector coordinates instead");
        }
    }
    /// Constructor with a 3D vector
    Axis(double a, double b, double c): vector_(vector3d(a,b,c)/norm(vector3d(a,b,c))) {}

    double projection(const Vector3D & positions) {
        auto dot_product = dot(vector_,positions);
        return abs(dot_product);
    }

private:
    /// axis coordinates
    Vector3D vector_;
};

#endif
