// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_DATA_STRUCTURES_C_OBJECT_ALLOCATION_SERVER_LOCALIZED_SEPARATOR_LINK_H_
#define SRC_C_INTERFACE_DATA_STRUCTURES_C_OBJECT_ALLOCATION_SERVER_LOCALIZED_SEPARATOR_LINK_H_

#include "src/data_structures/object_allocation_server.h"
#include "src/planar_reconstruction/localized_separator_link.h"

extern "C" {

struct c_ObjServer_LocSepLink {
  IRL::ObjectAllocationServer<IRL::LocalizedSeparatorLink>* obj_ptr = nullptr;
};

void c_ObjServer_LocSepLink_new(c_ObjServer_LocSepLink* a_self,
                                const std::size_t* a_number_to_allocate);

void c_ObjServer_LocSepLink_delete(c_ObjServer_LocSepLink* a_self);
}

#endif  // SRC_C_INTERFACE_DATA_STRUCTURES_C_OBJECT_ALLOCATION_SERVER_LOCALIZED_SEPARATOR_LINK_H_
