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
    enum Type {
        Linear,
        Radial,
    };

    /// Constructor with a 3D vector
    Axis(double a, double b, double c, Type type): vector_(vector3d(a, b, c)), type_(type) {
        /// normalize the axis
        if (vector_  == vector3d(0, 0, 0)) {
            throw CFilesError("Axis should not be null");
        } else {
            vector_ = vector_ / norm(vector_);
        }
    }

    static Axis parse(std::string string, Type type) {
        auto splitted = split(string,':');
        if (splitted.size() == 1) {
            if (splitted[0] == "X" or splitted[0] == "x") {
                return Axis(1, 0, 0, type);
            } else if (splitted[0] == "Y" or splitted[0] == "y") {
                return Axis(0, 1, 0, type);
            } else if (splitted[0] == "Z" or splitted[0] == "z") {
                return Axis(0, 0, 1, type);
            } else {
                throw CFilesError("Axis non implemented, enter vector coordinates instead");
            }
        } else if (splitted.size() == 3) {
            auto a = string2double(splitted[0]);
            auto b = string2double(splitted[1]);
            auto c = string2double(splitted[2]);
            return Axis(a, b, c, type);
        } else {
            throw CFilesError("Axis for density profile should be of size 3");
        }
    }

    /// Get the vector coordinates
    Vector3D& get_coordinates() { return vector_; }

    /// True if the axis type_ is 'Linear'
    bool is_linear() {return type_ == Linear;}

    /// projection on axis (may be negative)
    double projection(const Vector3D & positions) {
        return dot(vector_, positions);
    }

    /// radial distance to axis
    double radial(const Vector3D & positions) {
        auto dot_product = dot(vector_, positions);
        auto distance2 = norm2(positions) - dot_product * dot_product;
        return sqrt(distance2);
    }

private:
    /// axis coordinates
    Vector3D vector_;
    Type type_;
};

#endif
