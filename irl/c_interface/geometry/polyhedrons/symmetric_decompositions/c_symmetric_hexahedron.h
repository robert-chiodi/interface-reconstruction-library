// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_C_SYMMETRIC_HEXAHEDRON_H_
#define IRL_C_INTERFACE_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_C_SYMMETRIC_HEXAHEDRON_H_

#include "irl/geometry/polyhedrons/symmetric_decompositions/symmetric_hexahedron.h"

extern "C" {

struct c_SymHex {
  IRL::SymmetricHexahedron* obj_ptr = nullptr;
};

void c_SymHex_new(c_SymHex* a_self);
void c_SymHex_delete(c_SymHex* a_self);
void c_SymHex_construct(c_SymHex* a_self, const double* a_symmetric_hexahedron);
double c_SymHex_calculateVolume(const c_SymHex* a_self);
void c_SymHex_adjustCapToMatchVolume(c_SymHex* a_self,
                                     const double* a_correct_signed_volume);  
void c_SymHex_getBoundingPts(const c_SymHex* a_self, double* a_lower_pt,
                             double* a_upper_pt);
void c_SymHex_printToScreen(const c_SymHex* a_self);
void c_SymHex_getPt(const c_SymHex* a_self, const int* a_index, double* a_pt);
}

#endif // IRL_C_INTERFACE_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_C_SYMMETRIC_HEXAHEDRON_H_
