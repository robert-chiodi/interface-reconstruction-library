// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_STORED_PLANE_OF_EXISTENCE_TPP_
#define SRC_GEOMETRY_POLYGONS_STORED_PLANE_OF_EXISTENCE_TPP_

namespace IRL {

inline StoredPlaneOfExistence::StoredPlaneOfExistence(void)
    : plane_of_existence_m{Normal::fromScalarConstant(0.0), 0.0} {}

inline StoredPlaneOfExistence::StoredPlaneOfExistence(
    const Plane& a_plane_of_existence)
    : plane_of_existence_m(a_plane_of_existence) {}

inline void StoredPlaneOfExistence::setPlaneOfExistence_derived(
    const Plane& a_plane) {
  plane_of_existence_m = a_plane;
}

inline const Plane& StoredPlaneOfExistence::getPlaneOfExistence_derived(
    void) const {
  return plane_of_existence_m;
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYGONS_STORED_PLANE_OF_EXISTENCE_TPP_
