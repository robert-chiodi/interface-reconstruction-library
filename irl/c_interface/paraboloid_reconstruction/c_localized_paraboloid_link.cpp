// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/paraboloid_reconstruction/c_localized_paraboloid_link.h"

#include <cassert>

extern "C" {

void c_LocParabLink_new(c_LocParabLink* a_self, const c_PlanarLoc* a_localizer,
                        const c_Paraboloid* a_paraboloid) {
  assert(a_self->obj_ptr == nullptr);
  assert(a_localizer != nullptr);
  assert(a_localizer->obj_ptr != nullptr);
  assert(a_paraboloid != nullptr);
  assert(a_paraboloid->obj_ptr != nullptr);
  a_self->is_owning = true;
  a_self->obj_ptr = new IRL::LocalizedParaboloidLink<double>(
      a_localizer->obj_ptr, a_paraboloid->obj_ptr);
}

void c_LocParabLink_newFromObjectAllocationServer(
    c_LocParabLink* a_self,
    c_ObjServer_LocParabLink* a_object_allocation_server,
    const c_PlanarLoc* a_localizer, const c_Paraboloid* a_paraboloid) {
  assert(a_self->obj_ptr == nullptr);
  assert(a_object_allocation_server != nullptr);
  assert(a_object_allocation_server->obj_ptr != nullptr);
  assert(a_localizer != nullptr);
  assert(a_localizer->obj_ptr != nullptr);
  assert(a_paraboloid != nullptr);
  assert(a_paraboloid->obj_ptr != nullptr);
  a_self->is_owning = false;
  a_self->obj_ptr = a_object_allocation_server->obj_ptr->getNewObject();
  *a_self->obj_ptr = IRL::LocalizedParaboloidLink<double>(
      a_localizer->obj_ptr, a_paraboloid->obj_ptr);
}

void c_LocParabLink_delete(c_LocParabLink* a_self) {
  if (a_self->is_owning) {
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}

void c_LocParabLink_setId(c_LocParabLink* a_self, const int* a_id) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_id >= 0);
  a_self->obj_ptr->setId(static_cast<IRL::UnsignedIndex_t>(*a_id));
}

int c_LocParabLink_getId(const c_LocParabLink* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->getId());
}

void c_LocParabLink_setEdgeConnectivity(
    c_LocParabLink* a_self, const int* a_plane_index,
    const c_LocParabLink* a_ptr_to_neighbor) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_ptr_to_neighbor != nullptr);
  assert(a_ptr_to_neighbor->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  a_self->obj_ptr->setEdgeConnectivity(
      static_cast<IRL::UnsignedIndex_t>(*a_plane_index),
      a_ptr_to_neighbor->obj_ptr);
}

void c_LocParabLink_setEdgeConnectivityNull(c_LocParabLink* a_self,
                                            const int* a_plane_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_plane_index >= 0);
  a_self->obj_ptr->setEdgeConnectivity(
      static_cast<IRL::UnsignedIndex_t>(*a_plane_index), nullptr);
}

void c_LocParabLink_printToScreen(const c_LocParabLink* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  std::cout << *a_self->obj_ptr << std::endl;
}
}
