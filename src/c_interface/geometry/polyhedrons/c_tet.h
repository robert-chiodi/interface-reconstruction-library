// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_TET_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_TET_H_

#include "src/geometry/polyhedrons/tet.h"

extern "C" {

struct c_Tet {
  IRL::Tet* obj_ptr = nullptr;
};

void c_Tet_new(c_Tet* a_ptr);

void c_Tet_delete(c_Tet* a_ptr);

void c_Tet_construct(c_Tet* a_ptr, const double* a_tet_pts);

double c_Tet_calculateVolume(const c_Tet* a_self);

void c_Tet_getBoundingPts(c_Tet* a_ptr, double* a_lower_pt, double* a_upper_pt);

void c_Tet_printToScreen(const c_Tet* a_self);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_TET_H_
