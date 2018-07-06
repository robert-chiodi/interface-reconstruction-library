// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_C_SYMMETRIC_PYRAMID_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_C_SYMMETRIC_PYRAMID_H_

#include "src/geometry/polyhedrons/symmetric_decompositions/symmetric_pyramid.h"

extern "C" {

struct c_SymPyrmd {
  IRL::SymmetricPyramid* obj_ptr = nullptr;
};

void c_SymPyrmd_new(c_SymPyrmd* a_self);
void c_SymPyrmd_delete(c_SymPyrmd* a_self);
void c_SymPyrmd_construct(c_SymPyrmd* a_self, const double* a_symmetric_pyramid);
double c_SymPyrmd_calculateVolume(const c_SymPyrmd* a_self);
void c_SymPyrmd_getBoundingPts(const c_SymPyrmd* a_self, double* a_lower_pt,
                             double* a_upper_pt);
void c_SymPyrmd_printToScreen(const c_SymPyrmd* a_self);
void c_SymPyrmd_getPt(const c_SymPyrmd* a_self, const int* a_index, double* a_pt);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_C_SYMMETRIC_PYRAMID_H_
