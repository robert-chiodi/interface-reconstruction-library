// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/geometry/polyhedrons/c_dodecahedron.h"

#include <cassert>
#include <cstring>
#include <iostream>

extern "C" {

void c_Dod_new(c_Dod* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::Dodecahedron;
}

void c_Dod_delete(c_Dod* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_Dod_construct(c_Dod* a_self, const double* a_Dodecahedron) {
  assert(a_self != nullptr);
  (*a_self->obj_ptr) =
      IRL::Dodecahedron::fromRawDoublePointer(8, a_Dodecahedron);
}

double c_Dod_calculateVolume(const c_Dod* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return a_self->obj_ptr->calculateVolume();
}

void c_Dod_printToScreen(const c_Dod* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  std::cout << *a_self->obj_ptr << std::endl;
}

void c_Dod_getBoundingPts(c_Dod* a_self, double* a_lower_pt,
                          double* a_upper_pt) {
  assert(a_self != nullptr);
  const auto& lower_pt = a_self->obj_ptr->getLowerLimits();
  const auto& upper_pt = a_self->obj_ptr->getUpperLimits();
  std::memcpy(&a_lower_pt[0], &lower_pt[0], 3 * sizeof(double));
  std::memcpy(&a_upper_pt[0], &upper_pt[0], 3 * sizeof(double));
}
}
