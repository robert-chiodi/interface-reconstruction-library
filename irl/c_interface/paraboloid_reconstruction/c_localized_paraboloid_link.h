// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_PARABOLOID_RECONSTRUCTION_C_LOCALIZED_PARABOLOID_LINK_H_
#define IRL_C_INTERFACE_PARABOLOID_RECONSTRUCTION_C_LOCALIZED_PARABOLOID_LINK_H_

#include "irl/c_interface/data_structures/c_object_allocation_server_localized_paraboloid_link.h"
#include "irl/c_interface/paraboloid_reconstruction/c_paraboloid.h"
#include "irl/c_interface/planar_reconstruction/c_localizers.h"
#include "irl/data_structures/object_allocation_server.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

extern "C" {

struct c_LocParabLink {
  IRL::LocalizedParaboloidLink<double>* obj_ptr;
  bool is_owning = false;
};

void c_LocParabLink_new(c_LocParabLink* a_self, const c_PlanarLoc* a_localizer,
                        const c_Paraboloid* a_paraboloid);

void c_LocParabLink_newFromObjectAllocationServer(
    c_LocParabLink* a_self,
    c_ObjServer_LocParabLink* a_object_allocation_server,
    const c_PlanarLoc* a_localizer, const c_Paraboloid* a_paraboloid);

void c_LocParabLink_delete(c_LocParabLink* a_self);

void c_LocParabLink_setId(c_LocParabLink* a_self, const int* a_id);

int c_LocParabLink_getId(const c_LocParabLink* a_self);

void c_LocParabLink_setEdgeConnectivity(
    c_LocParabLink* a_self, const int* a_plane_index,
    const c_LocParabLink* a_ptr_to_neighbor);

void c_LocParabLink_setEdgeConnectivityNull(c_LocParabLink* a_self,
                                            const int* a_plane_index);

void c_LocParabLink_printToScreen(const c_LocParabLink* a_self);
}

#endif  // IRL_C_INTERFACE_PARABOLOID_RECONSTRUCTION_C_LOCALIZED_PARABOLOID_LINK_H_
