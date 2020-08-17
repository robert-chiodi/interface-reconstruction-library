// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/geometry/polyhedrons/c_pyramid.h"

#include <cassert>
#include <cstring>
#include <iostream>

extern "C" {

void c_Pyrmd_new(c_Pyrmd* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::Pyramid;
}

void c_Pyrmd_delete(c_Pyrmd* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_Pyrmd_construct(c_Pyrmd* a_self, const double* a_pyramid) {
  assert(a_self != nullptr);
  (*a_self->obj_ptr) =
      IRL::Pyramid::fromRawDoublePointer(5, a_pyramid);
}

double c_Pyrmd_calculateVolume(const c_Pyrmd* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return a_self->obj_ptr->calculateVolume();
}

void c_Pyrmd_printToScreen(const c_Pyrmd* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  std::cout << *a_self->obj_ptr << std::endl;
}

void c_Pyrmd_getBoundingPts(c_Pyrmd* a_self, double* a_lower_pt,
                          double* a_upper_pt) {
  assert(a_self != nullptr);
  const auto& lower_pt = a_self->obj_ptr->getLowerLimits();
  const auto& upper_pt = a_self->obj_ptr->getUpperLimits();
  std::memcpy(&a_lower_pt[0], &lower_pt[0], 3 * sizeof(double));
  std::memcpy(&a_upper_pt[0], &upper_pt[0], 3 * sizeof(double));
}
}
