// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_SEPARATORS_H_
#define SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_SEPARATORS_H_

#include "src/c_interface/data_structures/c_object_allocation_server_planar_separator.h"
#include "src/data_structures/object_allocation_server.h"
#include "src/geometry/general/normal.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/planar_reconstruction/planar_separator.h"

extern "C" {

struct c_PlanarSep {
  IRL::PlanarSeparator* obj_ptr = nullptr;
  bool is_owning = false;
};

void c_PlanarSep_new(c_PlanarSep* a_self);

void c_PlanarSep_newFromObjectAllocationServer(
    c_PlanarSep* a_self, c_ObjServer_PlanarSep* a_object_allocation_server);

void c_PlanarSep_delete(c_PlanarSep* a_self);

void c_PlanarSep_addPlane(c_PlanarSep* a_self, const double* a_normal,
                          const double* a_distance);

void c_PlanarSep_setNumberOfPlanes(c_PlanarSep* a_self,
                                   const int* a_number_to_set);

void c_PlanarSep_setPlane(c_PlanarSep* a_self, const int* a_plane_index_to_set,
                          const double* a_normal, const double* a_distance);

void c_PlanarSep_copy(c_PlanarSep* a_self,
                      const c_PlanarSep* a_other_planar_separator);

int c_PlanarSep_getNumberOfPlanes(const c_PlanarSep* a_self);

void c_PlanarSep_getPlane(c_PlanarSep* a_self, const int* a_index,
                          double* a_plane_listed);

bool c_PlanarSep_isFlipped(const c_PlanarSep* a_self);

void c_PlanarSep_printToScreen(const c_PlanarSep* a_self);

}  // end extern C

#endif  // SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_SEPARATORS_H_
