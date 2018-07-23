// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CFILES_AXIS_HPP
#define CFILES_AXIS_HPP

#include <cmath>
#include <sstream>

#include <chemfiles.hpp>

#include "Errors.hpp"
#include "utils.hpp"

using namespace chemfiles;

class Axis {
public:
    enum Type {
        Linear,
        Radial,
    };

    /// Constructor with a 3D vector
    Axis(double a, double b, double c, Type type): vector_(Vector3D(a, b, c)), type_(type) {
        /// normalize the axis
        if (vector_  == Vector3D(0, 0, 0)) {
            throw cfiles_error("Axis should not be null");
        } else {
            vector_ = vector_ / vector_.norm();
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
                throw cfiles_error("Unknown axis specification '{}'. It should be x, y, z or a:b:c", splitted[0]);
            }
        } else if (splitted.size() == 3) {
            auto a = string2double(splitted[0]);
            auto b = string2double(splitted[1]);
            auto c = string2double(splitted[2]);
            return Axis(a, b, c, type);
        } else {
            throw cfiles_error("Axis option should be x, y, z or a:b:c");
        }
    }

    /// Get a string describing the axis
    std::string str() const {
        if (vector_ == Vector3D(1, 0, 0)) {
            return "x";
        } else if (vector_ == Vector3D(0, 1, 0)) {
            return "y";
        } else if (vector_ == Vector3D(0, 0, 1)) {
            return "z";
        } else {
            std::stringstream ss;
            ss << "(" << vector_[0] << ", " << vector_[1] << ", " << vector_[2] << ")";
            return ss.str();
        }
    }

    /// Get the coordinates of the axis
    const Vector3D& vector() const { return vector_; }

    /// Check if the axis type is 'Linear'
    bool is_linear() const {return type_ == Linear;}

    /// Check if the axis type is 'Radial'
    bool is_radial() const {return type_ == Radial;}

    /// Project the given `point` on this axis. For radial axis, this returns
    /// the radial distance to the axis. For linear axis, the projection may be
    /// negative.
    double projection(const Vector3D& point) const {
        switch(type_) {
        case Linear:
            return dot(vector_, point);
        case Radial:
            return sqrt(point.norm() * point.norm() - dot(vector_, point) * dot(vector_, point));
        }
    }

private:
    Vector3D vector_;
    Type type_;
};

#endif
