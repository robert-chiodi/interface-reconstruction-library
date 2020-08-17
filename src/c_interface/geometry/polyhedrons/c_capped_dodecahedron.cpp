// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/geometry/polyhedrons/c_capped_dodecahedron.h"

extern "C" {

void c_CapDod_new(c_CapDod *a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::CappedDodecahedron;
}

void c_CapDod_delete(c_CapDod *a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_CapDod_construct(c_CapDod *a_self, const double *a_capped_dodecahedron) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  *a_self->obj_ptr =
      IRL::CappedDodecahedron::fromRawDoublePointer(9, a_capped_dodecahedron);
}

double c_CapDod_calculateVolume(const c_CapDod *a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return a_self->obj_ptr->calculateVolume();
}

void c_CapDod_calculateCentroid(c_CapDod *a_ptr, double *a_centroid) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  IRL::Pt centroid = a_ptr->obj_ptr->calculateCentroid();
  for (IRL::UnsignedIndex_t n = 0; n < 3; ++n) {
    a_centroid[n] = centroid[n];
  }
}

void c_CapDod_adjustCapToMatchVolume(c_CapDod *a_self,
                                     const double *a_correct_signed_volume) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->adjustCapToMatchVolume(*a_correct_signed_volume);
}

void c_CapDod_getBoundingPts(const c_CapDod *a_self, double *a_lower_pt,
                             double *a_upper_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  const IRL::Pt &lower_pt = a_self->obj_ptr->getLowerLimits();
  const IRL::Pt &upper_pt = a_self->obj_ptr->getUpperLimits();
  std::memcpy(&a_lower_pt[0], &lower_pt[0], 3 * sizeof(double));
  std::memcpy(&a_upper_pt[0], &upper_pt[0], 3 * sizeof(double));
}

void c_CapDod_getPt(const c_CapDod *a_self, const int *a_index, double *a_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0 &&
         *a_index < static_cast<int>(a_self->obj_ptr->getNumberOfVertices()));
  a_pt[0] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)].x();
  a_pt[1] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)].y();
  a_pt[2] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)].z();
}
}
