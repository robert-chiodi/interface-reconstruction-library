// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZED_SEPARATOR_GROUP_LINK_H_
#define SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZED_SEPARATOR_GROUP_LINK_H_

#include "src/c_interface/data_structures/c_object_allocation_server_localized_separator_link.h"
#include "src/c_interface/planar_reconstruction/c_localizers.h"
#include "src/c_interface/planar_reconstruction/c_planar_separator_path_group.h"
#include "src/data_structures/object_allocation_server.h"
#include "src/planar_reconstruction/localized_separator_group_link.h"

extern "C" {

struct c_LocSepGroupLink {
  IRL::LocalizedSeparatorGroupLink* obj_ptr;
  bool is_owning = false;
};

void c_LocSepGroupLink_new(c_LocSepGroupLink* a_self, const c_PlanarLoc* a_localizer,
                      const c_PlanarSepPathGroup* a_separator);

void c_LocSepGroupLink_delete(c_LocSepGroupLink* a_self);

void c_LocSepGroupLink_setId(c_LocSepGroupLink* a_self, const int* a_id);

int c_LocSepGroupLink_getId(const c_LocSepGroupLink* a_self);

void c_LocSepGroupLink_setEdgeConnectivity(c_LocSepGroupLink* a_self,
                                      const int* a_plane_index,
                                      const c_LocSepGroupLink* a_ptr_to_neighbor);

void c_LocSepGroupLink_setEdgeConnectivityNull(c_LocSepGroupLink* a_self,
                                          const int* a_plane_index);
}

#endif  // SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZED_SEPARATOR_GROUP_LINK_H_
