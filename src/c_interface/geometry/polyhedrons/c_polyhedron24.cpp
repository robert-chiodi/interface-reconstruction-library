// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/geometry/polyhedrons/c_polyhedron24.h"

extern "C" {

void c_Poly24_new(c_Poly24* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::Polyhedron24;
}

void c_Poly24_delete(c_Poly24* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_Poly24_construct(c_Poly24* a_self, const double* a_capped_dodecahedron) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr) =
      IRL::Polyhedron24::fromRawDoublePointer(14, a_capped_dodecahedron);
}

void c_Poly24_adjustCapToMatchVolume(c_Poly24* a_self,
                                     const double* a_correct_signed_volume) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->adjustCapToMatchVolume(*a_correct_signed_volume);
}

void c_Poly24_getBoundingPts(c_Poly24* a_self, double* a_lower_pt,
                             double* a_upper_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  IRL::Pt lower_pt = a_self->obj_ptr->getLowerLimits();
  IRL::Pt upper_pt = a_self->obj_ptr->getUpperLimits();
  std::memcpy(&a_lower_pt[0], &lower_pt[0], 3 * sizeof(double));
  std::memcpy(&a_upper_pt[0], &upper_pt[0], 3 * sizeof(double));
}

void c_Poly24_getPt(c_Poly24* a_self, const int* a_index, double* a_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0 &&
         *a_index < static_cast<int>(a_self->obj_ptr->getNumberOfVertices()));
  a_pt[0] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
                .getPt()
                .x();
  a_pt[1] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
                .getPt()
                .y();
  a_pt[2] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
                .getPt()
                .z();
}

void c_Poly24_setPt(c_Poly24* a_self, const int* a_index, const double* a_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0 &&
         *a_index < static_cast<int>(a_self->obj_ptr->getNumberOfVertices()));
  (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)][0] = a_pt[0];
  (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)][1] = a_pt[1];
  (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)][2] = a_pt[2];
}
}
