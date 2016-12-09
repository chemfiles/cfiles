// cfiles, an analysis frontend for the Chemfiles library
// Copyright (C) 2015-2016 Guillaume Fraux
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/

#ifndef CFILES_GEOMETRY_HPP
#define CFILES_GEOMETRY_HPP

#include <chemfiles.hpp>

namespace {
    using chemfiles::Vector3D;
    using chemfiles::dot;
    using chemfiles::cross;
    using chemfiles::norm;

    /// Compute the angle created by the vectors `r12` and `r23`.
    inline double angle(const Vector3D& r21, const Vector3D& r23) {
        double cos = dot(r21, r23) / (norm(r21) * norm(r23));
        cos = std::max(-1., std::min(1., cos));
        return std::acos(cos);
    }

    /// Compute the dihedral angle created by the vectors `r12`, `r23`, and `r34`.
    inline double dihedral(const Vector3D& r12, const Vector3D& r23, const Vector3D& r34) {
        auto a = cross(r12, r23);
        auto b = cross(r23, r34);
        return std::atan2(norm(r23) * dot(b, r12), dot(a, b));
    }
}


#endif
