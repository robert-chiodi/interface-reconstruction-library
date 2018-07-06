// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_OCTAHEDRON_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_OCTAHEDRON_H_

#include "src/geometry/polyhedrons/octahedron.h"

extern "C" {

struct c_Octa {
  IRL::Octahedron* obj_ptr = nullptr;
};

void c_Octa_new(c_Octa* a_self);
void c_Octa_delete(c_Octa* a_self);
void c_Octa_construct(c_Octa* a_self, const double* a_octahedron);
double c_Octa_calculateVolume(const c_Octa* a_self);
void c_Octa_printToScreen(const c_Octa* a_self);
void c_Octa_getBoundingPts(c_Octa* a_self, double* a_lower_pt,
                           double* a_upper_pt);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_OCTAHEDRON_H_
