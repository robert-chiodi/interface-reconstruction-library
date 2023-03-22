// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_PARABOLOID_RECONSTRUCTION_C_PARABOLOID_H_
#define IRL_C_INTERFACE_PARABOLOID_RECONSTRUCTION_C_PARABOLOID_H_

#include "irl/c_interface/data_structures/c_object_allocation_server_paraboloid.h"
#include "irl/c_interface/geometry/polygons/c_polygon.h"
#include "irl/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "irl/c_interface/paraboloid_reconstruction/c_triangulated_paraboloid.h"
#include "irl/c_interface/planar_reconstruction/c_separators.h"
#include "irl/data_structures/object_allocation_server.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/planar_reconstruction/planar_separator.h"

extern "C" {

struct c_Paraboloid {
  IRL::Paraboloid* obj_ptr = nullptr;
  bool is_owning = false;
};

void c_Paraboloid_new(c_Paraboloid* a_self);

void c_Paraboloid_newFromObjectAllocationServer(
    c_Paraboloid* a_self, c_ObjServer_Paraboloid* a_object_allocation_server);

void c_Paraboloid_delete(c_Paraboloid* a_self);

void c_Paraboloid_setDatum(c_Paraboloid* a_self, const double* a_datum);

void c_Paraboloid_setReferenceFrame(c_Paraboloid* a_self,
                                    const double* a_normal1,
                                    const double* a_normal2,
                                    const double* a_normal3);

void c_Paraboloid_setAlignedParaboloid(c_Paraboloid* a_self,
                                       const double* a_coeff_a,
                                       const double* a_coeff_b);

void c_Paraboloid_setParaboloidFromPolygon(c_Paraboloid* a_self,
                                           c_Poly* a_polygon);

void c_Paraboloid_setParaboloidFromPlanarSep(c_Paraboloid* a_self,
                                             c_PlanarSep* a_plane,
                                             double* a_pt_ref);

void c_Paraboloid_setParaboloidJibben(c_Paraboloid* a_self,
                                      c_PlanarSep* a_plane, c_RectCub* a_cell,
                                      int* a_npoly, double* a_vfrac,
                                      int* a_nvert, double* a_vert_coords);

void c_Paraboloid_setParaboloidFull(c_Paraboloid* a_self);

void c_Paraboloid_setParaboloidEmpty(c_Paraboloid* a_self);

void c_Paraboloid_copy(c_Paraboloid* a_self,
                       const c_Paraboloid* a_other_planar_separator);

void c_Paraboloid_getDatum(c_Paraboloid* a_self, double* a_datum);

void c_Paraboloid_getReferenceFrame(c_Paraboloid* a_self, double* a_frame);

void c_Paraboloid_getAlignedParaboloid(c_Paraboloid* a_self,
                                       double* a_aligned_paraboloid);

void c_Paraboloid_triangulateInsideCuboid(c_Paraboloid* a_self,
                                          c_RectCub* a_cell,
                                          c_TriangulatedParaboloid* a_surface);

double c_Paraboloid_getMeanCurvature(c_Paraboloid* a_self, c_RectCub* a_cell);

double c_Paraboloid_getSurfaceArea(c_Paraboloid* a_self, c_RectCub* a_cell);

void c_Paraboloid_printToScreen(const c_Paraboloid* a_self);

}  // end extern C

#endif  // IRL_C_INTERFACE_PARABOLOID_RECONSTRUCTION_C_PARABOLOID_H_
