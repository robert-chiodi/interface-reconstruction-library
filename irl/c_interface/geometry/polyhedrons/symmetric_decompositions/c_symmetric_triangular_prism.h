// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_C_SYMMETRIC_TRIANGULAR_PRISM_H_
#define IRL_C_INTERFACE_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_C_SYMMETRIC_TRIANGULAR_PRISM_H_

#include "irl/geometry/polyhedrons/symmetric_decompositions/symmetric_triangular_prism.h"

extern "C" {

struct c_SymTriPrism {
  IRL::SymmetricTriangularPrism* obj_ptr = nullptr;
};

void c_SymTriPrism_new(c_SymTriPrism* a_self);
void c_SymTriPrism_delete(c_SymTriPrism* a_self);
void c_SymTriPrism_construct(c_SymTriPrism* a_self, const double* a_symmetric_tri_prism);
double c_SymTriPrism_calculateVolume(const c_SymTriPrism* a_self);
void c_SymTriPrism_adjustCapToMatchVolume(c_SymTriPrism* a_self,
                                     const double* a_correct_signed_volume);  
void c_SymTriPrism_getBoundingPts(const c_SymTriPrism* a_self, double* a_lower_pt,
                             double* a_upper_pt);
void c_SymTriPrism_printToScreen(const c_SymTriPrism* a_self);
void c_SymTriPrism_getPt(const c_SymTriPrism* a_self, const int* a_index, double* a_pt);
}

#endif // IRL_C_INTERFACE_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_C_SYMMETRIC_TRIANGULAR_PRISM_H_
