// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_RECTANGULAR_CUBOID_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_RECTANGULAR_CUBOID_H_

#include "src/geometry/polyhedrons/rectangular_cuboid.h"

extern "C" {

struct c_RectCub {
  IRL::RectangularCuboid* obj_ptr = nullptr;
};

void c_RectCub_new(c_RectCub* a_self);

void c_RectCub_delete(c_RectCub* a_self);

void c_RectCub_construct(c_RectCub* a_self, const double* a_rectangular_cuboid);

void c_RectCub_construct_2pt(c_RectCub* a_self, const double* a_lower_pt,
                             const double* a_upper_pt);

double c_RectCub_calculateVolume(c_RectCub* a_self);

void c_RectCub_getBoundingPts(c_RectCub* a_self, double* a_lower_pt,
                              double* a_upper_pt);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_RECTANGULAR_CUBOID_H_
