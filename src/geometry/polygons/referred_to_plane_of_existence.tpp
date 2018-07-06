// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_REFERRED_TO_PLANE_OF_EXISTENCE_TPP_
#define SRC_GEOMETRY_POLYGONS_REFERRED_TO_PLANE_OF_EXISTENCE_TPP_

namespace IRL {

inline ReferredToPlaneOFExistence::ReferredToPlaneOFExistence(
    const Plane& a_plane_of_existence)
    : plane_of_existence_m(&a_plane_of_existence) {}

inline ReferredToPlaneOFExistence::ReferredToPlaneOFExistence(
    const Plane* a_plane_of_existence)
    : plane_of_existence_m(a_plane_of_existence) {}

inline const Plane& ReferredToPlaneOFExistence::getPlaneOfExistence_derived(
    void) const {
  assert(plane_of_existence_m != nullptr);
  return *plane_of_existence_m;
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYGONS_REFERRED_TO_PLANE_OF_EXISTENCE_TPP_
