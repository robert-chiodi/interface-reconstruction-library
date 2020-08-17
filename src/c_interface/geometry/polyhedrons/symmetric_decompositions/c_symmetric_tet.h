// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_C_SYMMETRIC_TET_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_C_SYMMETRIC_TET_H_

#include "src/geometry/polyhedrons/symmetric_decompositions/symmetric_tet.h"

extern "C" {

struct c_SymTet {
  IRL::SymmetricTet* obj_ptr = nullptr;
};

void c_SymTet_new(c_SymTet* a_self);
void c_SymTet_delete(c_SymTet* a_self);
void c_SymTet_construct(c_SymTet* a_self, const double* a_symmetric_tet);
double c_SymTet_calculateVolume(const c_SymTet* a_self);
void c_SymTet_getBoundingPts(const c_SymTet* a_self, double* a_lower_pt,
                             double* a_upper_pt);
void c_SymTet_printToScreen(const c_SymTet* a_self);
void c_SymTet_getPt(const c_SymTet* a_self, const int* a_index, double* a_pt);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_C_SYMMETRIC_TET_H_
