// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYGONS_C_POLYGON_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYGONS_C_POLYGON_H_

#include "src/c_interface/geometry/polygons/c_tri.h"
#include "src/c_interface/planar_reconstruction/c_localizers.h"

#include "src/geometry/polygons/polygon.h"

extern "C" {

struct c_Poly {
  IRL::Polygon* obj_ptr = nullptr;
};

void c_Poly_new(c_Poly* a_self);

void c_Poly_delete(c_Poly* a_self);

void c_Poly_construct(c_Poly* a_self, const int* a_number_of_pts,
                      const double* a_pts);

void c_Poly_calculateNormal(const c_Poly* a_self, double* a_normal);

void c_Poly_getLocalizer(const c_Poly* a_self, c_PlanarLoc* a_planar_localizer);

void c_Poly_reversePtOrdering(c_Poly* a_self);

void c_Poly_getBoundingPts(const c_Poly* a_self,
                           double* __restrict__ a_lower_pt,
                           double* __restrict__ a_upper_pt);

int c_Poly_getNumberOfPts(const c_Poly* a_self);

void c_Poly_getPt(const c_Poly* a_self, const int* a_index, double* a_pt);

int c_Poly_getNumberOfSimplicesInDecomposition(c_Poly* a_self);

void c_Poly_getSimplexFromDecomposition(c_Poly* a_self,
                                        const int* a_tri_number_to_get,
                                        c_Tri* a_tri_in_decomposition);

void c_Poly_zeroPolygon(c_Poly* a_self);

void c_Poly_calculateNearestPtOnSurface(const c_Poly* a_self,
                                        const double* a_pt,
                                        double* a_pt_on_polygon);

double c_Poly_calculateVolume(const c_Poly* a_self);

double c_Poly_calculateSign(const c_Poly* a_self);

void c_Poly_setPlaneOfExistence(c_Poly* a_self, const double* a_plane);

void c_Poly_calculateAndSetPlaneOfExistence(c_Poly* a_self);

void c_Poly_getPlaneOfExistence(const c_Poly* a_self, double* a_plane);

void c_Poly_calculateCentroid(const c_Poly* a_self, double* a_centroid);

void c_Poly_printToScreen(const c_Poly* a_self);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYGONS_C_POLYGON_H_
