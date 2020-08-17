// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_TRIANGULAR_PRISM_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_TRIANGULAR_PRISM_H_

#include "src/geometry/polyhedrons/triangular_prism.h"

extern "C" {

struct c_TriPrism {
  IRL::TriangularPrism* obj_ptr = nullptr;
};

void c_TriPrism_new(c_TriPrism* a_self);
void c_TriPrism_delete(c_TriPrism* a_self);
void c_TriPrism_construct(c_TriPrism* a_self, const double* a_triangular_prism);
double c_TriPrism_calculateVolume(const c_TriPrism* a_self);
void c_TriPrism_printToScreen(const c_TriPrism* a_self);
void c_TriPrism_getBoundingPts(c_TriPrism* a_self, double* a_lower_pt,
                          double* a_upper_pt);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_TRIANGULAR_PRISM_H_
