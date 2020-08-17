// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZERS_H_
#define SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZERS_H_

#include "src/c_interface/data_structures/c_object_allocation_server_planar_localizer.h"
#include "src/geometry/general/normal.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/planar_reconstruction/planar_localizer.h"

extern "C" {

struct c_PlanarLoc {
  IRL::PlanarLocalizer* obj_ptr = nullptr;
  bool is_owning = false;
};

void c_PlanarLoc_new(c_PlanarLoc* a_self);

void c_PlanarLoc_newFromObjectAllocationServer(
    c_PlanarLoc* a_self, c_ObjServer_PlanarLoc* a_object_allocation_server);

void c_PlanarLoc_delete(c_PlanarLoc* a_self);

void c_PlanarLoc_addPlane(c_PlanarLoc* a_self, const double* a_normal,
                          const double* a_distance);

void c_PlanarLoc_setNumberOfPlanes(c_PlanarLoc* a_self,
                                   const int* a_number_to_set);

void c_PlanarLoc_setPlane(c_PlanarLoc* a_self, const int* a_plane_index_to_set,
                          const double* a_normal, const double* a_distance);

void c_PlanarLoc_setFromRectangularCuboid(c_PlanarLoc* a_self,
                                          const double* a_lower_pt,
                                          const double* a_upper_pt);

void c_PlanarLoc_printToScreen(const c_PlanarLoc* a_self);
}
#endif  // SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZERS_H_
