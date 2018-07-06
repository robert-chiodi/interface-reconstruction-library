// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_DODECAHEDRON_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_DODECAHEDRON_H_

#include "src/geometry/polyhedrons/dodecahedron.h"

extern "C" {

struct c_Dod {
  IRL::Dodecahedron* obj_ptr = nullptr;
};

void c_Dod_new(c_Dod* a_self);
void c_Dod_delete(c_Dod* a_self);
void c_Dod_construct(c_Dod* a_self, const double* a_Dodecahedron);
double c_Dod_calculateVolume(const c_Dod* a_self);
void c_Dod_printToScreen(const c_Dod* a_self);
void c_Dod_getBoundingPts(c_Dod* a_self, double* a_lower_pt,
                          double* a_upper_pt);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_DODECAHEDRON_H_
