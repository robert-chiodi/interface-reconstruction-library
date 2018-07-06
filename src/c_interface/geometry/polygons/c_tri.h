// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYGONS_C_TRI_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYGONS_C_TRI_H_

#include "src/c_interface/planar_reconstruction/c_localizers.h"
#include "src/geometry/polygons/tri.h"

extern "C" {

struct c_Tri {
  IRL::Tri* obj_ptr = nullptr;
};

void c_Tri_new(c_Tri* a_ptr);

void c_Tri_delete(c_Tri* a_ptr);

void c_Tri_construct(c_Tri* a_ptr, const double* a_pts);

void c_Tri_getVertices(c_Tri* a_ptr, double* a_pts);

double c_Tri_calculateVolume(c_Tri* a_ptr);

void c_Tri_calculateCentroid(c_Tri* a_ptr, double* a_centroid);

void c_Tri_calculateNormal(c_Tri* a_ptr, double* a_normal);

void c_Tri_getLocalizer(const c_Tri* a_ptr, c_PlanarLoc* a_planar_localizer);

void c_Tri_reversePtOrdering(c_Tri* a_ptr);

void c_Tri_getBoundingPts(c_Tri* a_ptr, double* __restrict__ a_lower_pt,
                          double* __restrict__ a_upper_pt);

double c_Tri_calculateSign(const c_Tri* a_ptr);

void c_Tri_setPlaneOfExistence(c_Tri* a_ptr, const double* a_plane);

void c_Tri_calculateAndSetPlaneOfExistence(c_Tri* a_ptr);

void c_Tri_getPlaneOfExistence(const c_Tri* a_ptr, double* a_plane);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYGONS_C_TRI_H_
