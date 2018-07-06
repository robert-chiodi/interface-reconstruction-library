// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_PYRAMID_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_PYRAMID_H_

#include "src/geometry/polyhedrons/pyramid.h"

extern "C" {

struct c_Pyrmd {
  IRL::Pyramid* obj_ptr = nullptr;
};

void c_Pyrmd_new(c_Pyrmd* a_self);
void c_Pyrmd_delete(c_Pyrmd* a_self);
void c_Pyrmd_construct(c_Pyrmd* a_self, const double* a_octahedron);
double c_Pyrmd_calculateVolume(const c_Pyrmd* a_self);
void c_Pyrmd_printToScreen(const c_Pyrmd* a_self);
void c_Pyrmd_getBoundingPts(c_Pyrmd* a_self, double* a_lower_pt,
                          double* a_upper_pt);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_PYRAMID_H_
