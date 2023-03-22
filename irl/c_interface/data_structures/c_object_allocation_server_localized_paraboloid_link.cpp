// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/data_structures/c_object_allocation_server_localized_paraboloid_link.h"

#include <cassert>

extern "C" {

void c_ObjServer_LocParabLink_new(
    c_ObjServer_LocParabLink* a_self,
    const IRL::LargeOffsetIndex_t* a_number_to_allocate) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr =
      new IRL::ObjectAllocationServer<IRL::LocalizedParaboloidLink<double>>(
          *a_number_to_allocate);
}

void c_ObjServer_LocParabLink_delete(c_ObjServer_LocParabLink* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}
}
