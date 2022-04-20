// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <iostream>

#include "examples/simple_paraboloid_advector/basic_mesh.h"

void BasicMesh::setCellBoundaries(const IRL::Pt& a_bottom_bounding_box,
                                  const IRL::Pt& a_top_bounding_box) {
  // X direction
  dx_m = (a_top_bounding_box.x() - a_bottom_bounding_box.x()) /
         static_cast<double>(this->getNx());
  for (int i = -this->getNgc(); i < this->getNx() + 1 + this->getNgc(); ++i) {
    this->x(i) = a_bottom_bounding_box.x() + static_cast<double>(i) * dx_m;
  }

  // Y direction
  dy_m = (a_top_bounding_box.y() - a_bottom_bounding_box.y()) /
         static_cast<double>(this->getNy());
  for (int j = -this->getNgc(); j < this->getNy() + 1 + this->getNgc(); ++j) {
    this->y(j) = a_bottom_bounding_box.y() + static_cast<double>(j) * dy_m;
  }

  // Z direction
  dz_m = (a_top_bounding_box.z() - a_bottom_bounding_box.z()) /
         static_cast<double>(this->getNz());
  for (int k = -this->getNgc(); k < this->getNz() + 1 + this->getNgc(); ++k) {
    this->z(k) = a_bottom_bounding_box.z() + static_cast<double>(k) * dz_m;
  }
}

void BasicMesh::getIndices(const IRL::Pt& a_location, int* a_indices) const {
  // Localize x index
  int i =
      static_cast<int>((a_location[0] - this->x(this->imin())) / this->dx()) +
      this->imin();
  int j =
      static_cast<int>((a_location[1] - this->y(this->jmin())) / this->dy()) +
      this->jmin();
  int k =
      static_cast<int>((a_location[2] - this->z(this->kmin())) / this->dz()) +
      this->kmin();
  a_indices[0] = i;
  a_indices[1] = j;
  a_indices[2] = k;
}
