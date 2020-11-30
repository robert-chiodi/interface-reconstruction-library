// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_POLYGONS_REFERRED_TO_PLANE_OF_EXISTENCE_H_
#define IRL_GEOMETRY_POLYGONS_REFERRED_TO_PLANE_OF_EXISTENCE_H_

#include "irl/geometry/general/plane.h"

namespace IRL {

class ReferredToPlaneOFExistence {
 public:
  ReferredToPlaneOFExistence(void) = delete;

  explicit ReferredToPlaneOFExistence(const Plane& a_plane_of_existence);

  explicit ReferredToPlaneOFExistence(const Plane* a_plane_of_existence);

  const Plane& getPlaneOfExistence_derived(void) const;

 private:
  const Plane* plane_of_existence_m;
};
}  // namespace IRL

#include "irl/geometry/polygons/referred_to_plane_of_existence.tpp"

#endif // IRL_GEOMETRY_POLYGONS_REFERRED_TO_PLANE_OF_EXISTENCE_H_
