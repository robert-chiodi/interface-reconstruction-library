// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZER_LINK_H_
#define SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZER_LINK_H_

#include "src/c_interface/data_structures/c_object_allocation_server_localizer_link.h"
#include "src/c_interface/planar_reconstruction/c_localizers.h"
#include "src/data_structures/object_allocation_server.h"
#include "src/planar_reconstruction/localizer_link.h"

extern "C" {

struct c_LocLink {
  IRL::LocalizerLink* obj_ptr = nullptr;
  bool is_owning = false;
};

void c_LocLink_new(c_LocLink* a_self, const c_PlanarLoc* a_localizer);

void c_LocLink_newFromObjectAllocationServer(
    c_LocLink* a_self, c_ObjServer_LocLink* a_object_allocation_server,
    const c_PlanarLoc* a_localizer);

void c_LocLink_delete(c_LocLink* a_self);

void c_LocLink_setId(c_LocLink* a_self, const int* a_id);

int c_LocLink_getId(const c_LocLink* a_self);

void c_LocLink_setEdgeConnectivity(c_LocLink* a_self, const int* a_plane_index,
                                   const c_LocLink* a_self_to_neighbor);

void c_LocLink_setEdgeConnectivityNull(c_LocLink* a_self,
                                       const int* a_plane_index);
}

#endif  // SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZER_LINK_H_
