// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_CAPPED_Octahedron_VARIATIONS_C_CAPPED_Octahedron_LLL_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_CAPPED_Octahedron_VARIATIONS_C_CAPPED_Octahedron_LLL_H_

#include "src/geometry/polyhedrons/capped_octahedron_variations/capped_octahedron_LLL.h"

extern "C" {

struct c_CapOcta_LLL {
  IRL::CappedOctahedron_LLL* obj_ptr = nullptr;
};

void c_CapOcta_LLL_new(c_CapOcta_LLL* a_self);
void c_CapOcta_LLL_delete(c_CapOcta_LLL* a_self);
void c_CapOcta_LLL_construct(c_CapOcta_LLL* a_self, const double* a_octahedron);
void c_CapOcta_LLL_adjustCapToMatchVolume(c_CapOcta_LLL* a_self,
                                     const double* a_correct_signed_volume);
double c_CapOcta_LLL_calculateVolume(const c_CapOcta_LLL* a_self);
void c_CapOcta_LLL_getBoundingPts(const c_CapOcta_LLL* a_self, double* a_lower_pt,
                             double* a_upper_pt);
void c_CapOcta_LLL_printToScreen(const c_CapOcta_LLL* a_self);
void c_CapOcta_LLL_getPt(const c_CapOcta_LLL* a_self, const int* a_index, double* a_pt);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_CAPPED_Octahedron_VARIATIONS_C_CAPPED_Octahedron_LLL_H_
