// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_POLYHEDRON24_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_POLYHEDRON24_H_

#include "src/geometry/polyhedrons/polyhedron_24.h"

extern "C" {

struct c_Poly24 {
  IRL::Polyhedron24* obj_ptr = nullptr;
};

void c_Poly24_new(c_Poly24* a_self);

void c_Poly24_delete(c_Poly24* a_self);

void c_Poly24_construct(c_Poly24* a_self, const double* a_dodecahedron);

void c_Poly24_adjustCapToMatchVolume(c_Poly24* a_self,
                                     const double* a_correct_signed_volume);

void c_Poly24_getBoundingPts(c_Poly24* a_self, double* a_lower_pt,
                             double* a_upper_pt);

void c_Poly24_getPt(c_Poly24* a_self, const int* a_index, double* a_pt);

void c_Poly24_setPt(c_Poly24* a_self, const int* a_index, const double* a_pt);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_POLYHEDRON24_H_
