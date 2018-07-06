// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/planar_reconstruction/c_localized_separator.h"

#include <cassert>

extern "C" {

void c_LocSep_new(c_LocSep* a_self, const c_PlanarLoc* a_localizer,
                  const c_PlanarSep* a_separator) {
  assert(a_self->obj_ptr == nullptr);
  assert(a_localizer != nullptr);
  assert(a_localizer->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  a_self->is_owning = true;
  a_self->obj_ptr =
      new IRL::LocalizedSeparator(a_localizer->obj_ptr, a_separator->obj_ptr);
}

void c_LocSep_newFromObjectAllocationServer(
    c_LocSep* a_self, c_ObjServer_LocSep* a_object_allocation_server,
    const c_PlanarLoc* a_localizer, const c_PlanarSep* a_separator) {
  assert(a_self->obj_ptr == nullptr);
  assert(a_object_allocation_server != nullptr);
  assert(a_object_allocation_server->obj_ptr != nullptr);
  assert(a_localizer != nullptr);
  assert(a_localizer->obj_ptr != nullptr);
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  a_self->is_owning = false;
  a_self->obj_ptr = a_object_allocation_server->obj_ptr->getNewObject();
  *a_self->obj_ptr =
      IRL::LocalizedSeparator(a_localizer->obj_ptr, a_separator->obj_ptr);
}

void c_LocSep_delete(c_LocSep* a_self) {
  if (a_self->is_owning) {
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}
}
