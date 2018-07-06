// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/geometry/polyhedrons/symmetric_decompositions/c_symmetric_triangular_prism.h"

extern "C" {

void c_SymTriPrism_new(c_SymTriPrism* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::SymmetricTriangularPrism;
}

void c_SymTriPrism_delete(c_SymTriPrism* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_SymTriPrism_construct(c_SymTriPrism* a_self, const double* a_symmetric_tri_prism) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  *a_self->obj_ptr =
      IRL::SymmetricTriangularPrism::fromRawDoublePointer(11, a_symmetric_tri_prism);
}

double c_SymTriPrism_calculateVolume(const c_SymTriPrism* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return a_self->obj_ptr->calculateVolume();
}

void c_SymTriPrism_adjustCapToMatchVolume(c_SymTriPrism* a_self,
                                     const double* a_correct_signed_volume) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->adjustCapToMatchVolume(*a_correct_signed_volume);
}  

void c_SymTriPrism_getBoundingPts(const c_SymTriPrism* a_self, double* a_lower_pt,
                             double* a_upper_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  const IRL::Pt& lower_pt = a_self->obj_ptr->getLowerLimits();
  const IRL::Pt& upper_pt = a_self->obj_ptr->getUpperLimits();
  std::memcpy(&a_lower_pt[0], &lower_pt[0], 3 * sizeof(double));
  std::memcpy(&a_upper_pt[0], &upper_pt[0], 3 * sizeof(double));
}

void c_SymTriPrism_printToScreen(const c_SymTriPrism* a_self){
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  std::cout << *a_self->obj_ptr <<std::endl;
}

void c_SymTriPrism_getPt(const c_SymTriPrism* a_self, const int* a_index, double* a_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0 &&
         *a_index < static_cast<int>(a_self->obj_ptr->getNumberOfVertices()));
  a_pt[0] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)].x();
  a_pt[1] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)].y();
  a_pt[2] = (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)].z();
}
}
