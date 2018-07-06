// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_STORED_PLANE_OF_EXISTENCE_H_
#define SRC_GEOMETRY_POLYGONS_STORED_PLANE_OF_EXISTENCE_H_

#include "src/geometry/general/plane.h"

namespace IRL {

class StoredPlaneOfExistence {
 public:
  StoredPlaneOfExistence(void);

  explicit StoredPlaneOfExistence(const Plane& a_plane_of_existence);
  void setPlaneOfExistence_derived(const Plane& a_plane);

  const Plane& getPlaneOfExistence_derived(void) const;

 private:
  Plane plane_of_existence_m;
};

}  // namespace IRL

#include "src/geometry/polygons/stored_plane_of_existence.tpp"

#endif  // SRC_GEOMETRY_POLYGONS_STORED_PLANE_OF_EXISTENCE_H_
