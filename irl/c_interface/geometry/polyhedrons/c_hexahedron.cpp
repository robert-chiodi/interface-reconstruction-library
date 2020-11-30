// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/geometry/polyhedrons/c_hexahedron.h"

extern "C" {

void c_Hex_new(c_Hex* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::Hexahedron;
}

void c_Hex_delete(c_Hex* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_Hex_construct(c_Hex* a_self,
                     const double* a_hexahedron) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr) =
      IRL::Hexahedron::fromRawDoublePointer(8, a_hexahedron);
}

double c_Hex_calculateVolume(const c_Hex* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return a_self->obj_ptr->calculateVolume();
}

void c_Hex_getVertices(const c_Hex* a_self, double* a_pts) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  for(IRL::UnsignedIndex_t v = 0; v < a_self->obj_ptr->getNumberOfVertices();++v){
    std::memcpy(&a_pts[3*v], &(*a_self->obj_ptr)[v][0], 3 * sizeof(double));
  }
}

void c_Hex_getBoundingPts(const c_Hex* a_self, double* a_lower_pt,
                              double* a_upper_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  IRL::Pt lower_pt = a_self->obj_ptr->getLowerLimits();
  IRL::Pt upper_pt = a_self->obj_ptr->getUpperLimits();
  std::memcpy(&a_lower_pt[0], &lower_pt[0], 3 * sizeof(double));
  std::memcpy(&a_upper_pt[0], &upper_pt[0], 3 * sizeof(double));
}

void c_Hex_printToScreen(const c_Hex* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  std::cout << *a_self->obj_ptr << std::endl;
}
  
}
