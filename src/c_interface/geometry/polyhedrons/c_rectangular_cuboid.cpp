// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"

extern "C" {

void c_RectCub_new(c_RectCub* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::RectangularCuboid;
}

void c_RectCub_delete(c_RectCub* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_RectCub_construct(c_RectCub* a_self,
                         const double* a_rectangular_cuboid) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr) =
      IRL::RectangularCuboid::fromRawDoublePointer(8, a_rectangular_cuboid);
}

void c_RectCub_construct_2pt(c_RectCub* a_self, const double* a_lower_pt,
                             const double* a_upper_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr) = IRL::RectangularCuboid::fromBoundingPts(
      IRL::Pt::fromRawDoublePointer(a_lower_pt),
      IRL::Pt::fromRawDoublePointer(a_upper_pt));
}

double c_RectCub_calculateVolume(c_RectCub* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return a_self->obj_ptr->calculateVolume();
}

void c_RectCub_getBoundingPts(c_RectCub* a_self, double* a_lower_pt,
                              double* a_upper_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  IRL::Pt lower_pt = a_self->obj_ptr->getLowerLimits();
  IRL::Pt upper_pt = a_self->obj_ptr->getUpperLimits();
  std::memcpy(&a_lower_pt[0], &lower_pt[0], 3 * sizeof(double));
  std::memcpy(&a_upper_pt[0], &upper_pt[0], 3 * sizeof(double));
}
}
