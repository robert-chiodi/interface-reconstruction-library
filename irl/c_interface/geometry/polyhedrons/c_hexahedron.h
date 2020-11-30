// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_HEXAHEDRON_H_
#define IRL_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_HEXAHEDRON_H_

#include "irl/geometry/polyhedrons/hexahedron.h"

extern "C" {

struct c_Hex {
  IRL::Hexahedron* obj_ptr = nullptr;
};

void c_Hex_new(c_Hex* a_self);

void c_Hex_delete(c_Hex* a_self);

void c_Hex_construct(c_Hex* a_self, const double* a_rectangular_cuboid);

double c_Hex_calculateVolume(const c_Hex* a_self);

void c_Hex_getVertices(const c_Hex* a_self, double* a_pts);

void c_Hex_getBoundingPts(const c_Hex* a_self, double* a_lower_pt,
                              double* a_upper_pt);

void c_Hex_printToScreen(const c_Hex* a_self);  
}

#endif // IRL_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_HEXAHEDRON_H_
