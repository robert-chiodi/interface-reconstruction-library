// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_PARABOLOID_RECONSTRUCTION_C_TRIANGULATED_PARABOLOID_H_
#define IRL_C_INTERFACE_PARABOLOID_RECONSTRUCTION_C_TRIANGULATED_PARABOLOID_H_

#include "irl/c_interface/data_structures/c_object_allocation_server_triangulated_paraboloid.h"
#include "irl/data_structures/object_allocation_server.h"
#include "irl/surface_mesher/triangulated_surface.h"

extern "C" {

struct c_TriangulatedParaboloid {
  IRL::TriangulatedSurfaceOutput* obj_ptr = nullptr;
  bool is_owning = false;
};

void c_TriangulatedParaboloid_new(c_TriangulatedParaboloid* a_self);

void c_TriangulatedParaboloid_newFromObjectAllocationServer(
    c_TriangulatedParaboloid* a_self,
    c_ObjServer_TriangulatedParaboloid* a_object_allocation_server);

void c_TriangulatedParaboloid_delete(c_TriangulatedParaboloid* a_self);

void c_TriangulatedParaboloid_clear(c_TriangulatedParaboloid* a_self);

int c_TriangulatedParaboloid_getNumberOfTriangles(
    const c_TriangulatedParaboloid* a_self);

void c_TriangulatedParaboloid_getPt(c_TriangulatedParaboloid* a_self,
                                    const int* a_tri_index,
                                    const int* a_vert_index, double* a_pt);

}  // end extern C

#endif  // IRL_C_INTERFACE_PARABOLOID_RECONSTRUCTION_C_TRIANGULATED_PARABOLOID_H_
