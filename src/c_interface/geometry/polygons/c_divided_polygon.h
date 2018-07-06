// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GEOMETRY_POLYGONS_C_DIVIDED_POLYGON_H_
#define SRC_C_INTERFACE_GEOMETRY_POLYGONS_C_DIVIDED_POLYGON_H_

#include "src/c_interface/geometry/polygons/c_polygon.h"
#include "src/c_interface/geometry/polygons/c_tri.h"
#include "src/c_interface/planar_reconstruction/c_localizers.h"
#include "src/geometry/polygons/divided_polygon.h"
#include "src/geometry/polygons/polygon.h"

extern "C" {

struct c_DivPoly {
  IRL::DividedPolygon* obj_ptr = nullptr;
};

void c_DivPoly_new(c_DivPoly* a_self);
void c_DivPoly_delete(c_DivPoly* a_self);
void c_DivPoly_construct(c_DivPoly* a_self, const int* a_number_of_pts,
                         const double* a_pts);
void c_DivPoly_constructFromPolygon(c_DivPoly* a_self, c_Poly* a_polygon_ptr);
void c_DivPoly_resetCentroid(c_DivPoly* a_self);
int c_DivPoly_getNumberOfSimplicesInDecomposition(c_DivPoly* a_self);
void c_DivPoly_getSimplexFromDecomposition(c_DivPoly* a_self,
                                           const int* a_tri_number_to_get,
                                           c_Tri* a_tri_in_decomposition);
void c_DivPoly_getLocalizer(const c_DivPoly* a_self,
                            c_PlanarLoc* a_planar_localizer);
void c_DivPoly_calculateNormal(c_DivPoly* a_self, double* a_normal);
void c_DivPoly_reversePtOrdering(c_DivPoly* a_self);
void c_DivPoly_getBoundingPts(c_DivPoly* a_self,
                              double* __restrict__ a_lower_pt,
                              double* __restrict__ a_upper_pt);
int c_DivPoly_getNumberOfPts(const c_DivPoly* a_self);
void c_DivPoly_getPt(const c_DivPoly* a_self, const int* a_index, double* a_pt);
void c_DivPoly_zeroPolygon(c_DivPoly* a_self);
double c_DivPoly_calculateSurfaceArea(const c_DivPoly* a_self);
double c_DivPoly_calculateSign(const c_DivPoly* a_self);
void c_DivPoly_setPlaneOfExistence(c_DivPoly* a_self, const double* a_plane);
void c_DivPoly_calculateAndSetPlaneOfExistence(c_DivPoly* a_self);
void c_DivPoly_getPlaneOfExistence(const c_DivPoly* a_self, double* a_plane);
void c_DivPoly_printToScreen(const c_DivPoly* a_self);
}

#endif  // SRC_C_INTERFACE_GEOMETRY_POLYGONS_C_DIVIDED_POLYGON_H_
