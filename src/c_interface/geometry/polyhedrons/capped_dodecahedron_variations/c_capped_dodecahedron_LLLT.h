// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_CAPPED_DODECAHEDRON_VARIATIONS_C_CAPPED_DODECAHEDRON_LLLT_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_CAPPED_DODECAHEDRON_VARIATIONS_C_CAPPED_DODECAHEDRON_LLLT_H_

#include "src/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_LLLT.h"

extern "C" {

struct c_CapDod_LLLT {
  IRL::CappedDodecahedron_LLLT* obj_ptr = nullptr;
};

void c_CapDod_LLLT_new(c_CapDod_LLLT* a_self);
void c_CapDod_LLLT_delete(c_CapDod_LLLT* a_self);
void c_CapDod_LLLT_construct(c_CapDod_LLLT* a_self, const double* a_dodecahedron);
void c_CapDod_LLLT_adjustCapToMatchVolume(c_CapDod_LLLT* a_self,
                                     const double* a_correct_signed_volume);
void c_CapDod_LLLT_getBoundingPts(const c_CapDod_LLLT* a_self, double* a_lower_pt,
                             double* a_upper_pt);
void c_CapDod_LLLT_printToScreen(const c_CapDod_LLLT* a_self);
void c_CapDod_LLLT_getPt(const c_CapDod_LLLT* a_self, const int* a_index, double* a_pt);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_CAPPED_DODECAHEDRON_VARIATIONS_C_CAPPED_DODECAHEDRON_LLLT_H_
