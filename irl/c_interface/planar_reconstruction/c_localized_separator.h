// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZED_SEPARATOR_H_
#define IRL_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZED_SEPARATOR_H_

#include "irl/c_interface/data_structures/c_object_allocation_server_localized_separator.h"
#include "irl/c_interface/planar_reconstruction/c_localizers.h"
#include "irl/c_interface/planar_reconstruction/c_separators.h"
#include "irl/data_structures/object_allocation_server.h"
#include "irl/planar_reconstruction/localized_separator.h"

extern "C" {

struct c_LocSep {
  IRL::LocalizedSeparator* obj_ptr;
  bool is_owning = false;
};

void c_LocSep_new(c_LocSep* a_self, const c_PlanarLoc* a_localizer,
                  const c_PlanarSep* a_separator);

void c_LocSep_newFromObjectAllocationServer(
    c_LocSep* a_self, c_ObjServer_LocSep* a_object_allocation_server,
    const c_PlanarLoc* a_localizer, const c_PlanarSep* a_separator);

void c_LocSep_delete(c_LocSep* a_self);
}

#endif // IRL_C_INTERFACE_PLANAR_RECONSTRUCTION_C_LOCALIZED_SEPARATOR_H_
