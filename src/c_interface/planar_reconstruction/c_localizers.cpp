// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/planar_reconstruction/c_localizers.h"

#include <cassert>
#include <iostream>

extern "C" {

void c_PlanarLoc_new(c_PlanarLoc* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->is_owning = true;
  a_self->obj_ptr = new IRL::PlanarLocalizer;
}

void c_PlanarLoc_newFromObjectAllocationServer(
    c_PlanarLoc* a_self, c_ObjServer_PlanarLoc* a_object_allocation_server) {
  assert(a_self->obj_ptr == nullptr);
  assert(a_object_allocation_server != nullptr);
  assert(a_object_allocation_server->obj_ptr != nullptr);
  a_self->is_owning = false;
  a_self->obj_ptr = a_object_allocation_server->obj_ptr->getNewObject();
}

void c_PlanarLoc_delete(c_PlanarLoc* a_self) {
  if (a_self->is_owning) {
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}

void c_PlanarLoc_addPlane(c_PlanarLoc* a_self, const double* a_normal,
                          const double* a_distance) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->addPlane(
      IRL::Plane(IRL::Normal::fromRawDoublePointer(a_normal), *a_distance));
}

void c_PlanarLoc_setNumberOfPlanes(c_PlanarLoc* a_self,
                                   const int* a_number_to_set) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_number_to_set >= 0);
  a_self->obj_ptr->setNumberOfPlanes(
      static_cast<IRL::UnsignedIndex_t>(*a_number_to_set));
}

void c_PlanarLoc_setPlane(c_PlanarLoc* a_self, const int* a_plane_index_to_set,
                          const double* a_normal, const double* a_distance) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_plane_index_to_set >= 0);
  (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_plane_index_to_set)] =
      IRL::Plane(IRL::Normal::fromRawDoublePointer(a_normal), *a_distance);
}

void c_PlanarLoc_setFromRectangularCuboid(c_PlanarLoc* a_self,
                                          const double* a_lower_pt,
                                          const double* a_upper_pt) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr) = IRL::RectangularCuboid::fromBoundingPts(
                           IRL::Pt::fromRawDoublePointer(a_lower_pt),
                           IRL::Pt::fromRawDoublePointer(a_upper_pt))
                           .getLocalizer();
}

void c_PlanarLoc_printToScreen(const c_PlanarLoc* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  std::cout << (*a_self->obj_ptr);
}
}
