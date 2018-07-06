// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZED_SEPARATOR_LINK_H_
#define SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZED_SEPARATOR_LINK_H_

#include "src/c_interface/data_structures/c_object_allocation_server_localized_separator_link.h"
#include "src/c_interface/planar_reconstruction/c_localizers.h"
#include "src/c_interface/planar_reconstruction/c_separators.h"
#include "src/data_structures/object_allocation_server.h"
#include "src/planar_reconstruction/localized_separator_link.h"

extern "C" {

struct c_LocSepLink {
  IRL::LocalizedSeparatorLink* obj_ptr;
  bool is_owning = false;
};

void c_LocSepLink_new(c_LocSepLink* a_self, const c_PlanarLoc* a_localizer,
                      const c_PlanarSep* a_separator);

void c_LocSepLink_newFromObjectAllocationServer(
    c_LocSepLink* a_self, c_ObjServer_LocSepLink* a_object_allocation_server,
    const c_PlanarLoc* a_localizer, const c_PlanarSep* a_separator);

void c_LocSepLink_delete(c_LocSepLink* a_self);

void c_LocSepLink_setId(c_LocSepLink* a_self, const int* a_id);

int c_LocSepLink_getId(const c_LocSepLink* a_self);

void c_LocSepLink_setEdgeConnectivity(c_LocSepLink* a_self,
                                      const int* a_plane_index,
                                      const c_LocSepLink* a_ptr_to_neighbor);

void c_LocSepLink_setEdgeConnectivityNull(c_LocSepLink* a_self,
                                          const int* a_plane_index);

void c_LocSepLink_printToScreen(const c_LocSepLink* a_self);
}

#endif  // SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZED_SEPARATOR_LINK_H_
