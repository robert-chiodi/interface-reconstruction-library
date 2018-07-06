// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/geometry/polyhedrons/c_tet.h"

extern "C" {

void c_Tet_new(c_Tet* a_ptr) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr == nullptr);
  a_ptr->obj_ptr = new IRL::Tet;
}

void c_Tet_delete(c_Tet* a_ptr) {
  delete a_ptr->obj_ptr;
  a_ptr->obj_ptr = nullptr;
}

void c_Tet_construct(c_Tet* a_ptr, const double* a_tet_pts) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  (*a_ptr->obj_ptr) = IRL::Tet::fromRawDoublePointer(4, a_tet_pts);
}

double c_Tet_calculateVolume(const c_Tet* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return a_self->obj_ptr->calculateVolume();
}

void c_Tet_getBoundingPts(c_Tet* a_ptr, double* a_lower_pt,
                          double* a_upper_pt) {
  assert(a_ptr != nullptr);
  assert(a_ptr->obj_ptr != nullptr);
  IRL::Pt lower_pt = a_ptr->obj_ptr->getLowerLimits();
  IRL::Pt upper_pt = a_ptr->obj_ptr->getUpperLimits();
  std::memcpy(&a_lower_pt[0], &lower_pt[0], 3 * sizeof(double));
  std::memcpy(&a_upper_pt[0], &upper_pt[0], 3 * sizeof(double));
}

void c_Tet_printToScreen(const c_Tet* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  std::cout << *a_self->obj_ptr << std::endl;
}

}
