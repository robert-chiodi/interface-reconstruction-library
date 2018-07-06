// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/planar_reconstruction/c_separators.h"

#include <iostream>

#include "src/interface_reconstruction_methods/volume_fraction_matching.h"
#include "src/parameters/constants.h"

extern "C" {

void c_PlanarSep_new(c_PlanarSep* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->is_owning = true;
  a_self->obj_ptr = new IRL::PlanarSeparator;
}

void c_PlanarSep_newFromObjectAllocationServer(
    c_PlanarSep* a_self, c_ObjServer_PlanarSep* a_object_allocation_server) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  assert(a_object_allocation_server != nullptr);
  assert(a_object_allocation_server->obj_ptr != nullptr);
  a_self->is_owning = false;
  a_self->obj_ptr = a_object_allocation_server->obj_ptr->getNewObject();
}

void c_PlanarSep_delete(c_PlanarSep* a_self) {
  if (a_self->is_owning) {
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}

void c_PlanarSep_addPlane(c_PlanarSep* a_self, const double* a_normal,
                          const double* a_distance) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->addPlane(
      IRL::Plane(IRL::Normal::fromRawDoublePointer(a_normal), *a_distance));
}

void c_PlanarSep_setNumberOfPlanes(c_PlanarSep* a_self,
                                   const int* a_number_to_set) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_number_to_set >= 0);
  a_self->obj_ptr->setNumberOfPlanes(
      static_cast<IRL::UnsignedIndex_t>(*a_number_to_set));
}

void c_PlanarSep_setPlane(c_PlanarSep* a_self, const int* a_plane_index_to_set,
                          const double* a_normal, const double* a_distance) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_plane_index_to_set >= 0);
  (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_plane_index_to_set)] =
      IRL::Plane(IRL::Normal::fromRawDoublePointer(a_normal), *a_distance);
}

void c_PlanarSep_copy(c_PlanarSep* a_self,
                      const c_PlanarSep* a_other_planar_separator) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_other_planar_separator != nullptr);
  assert(a_other_planar_separator->obj_ptr != nullptr);
  (*a_self->obj_ptr) = (*a_other_planar_separator->obj_ptr);
}

int c_PlanarSep_getNumberOfPlanes(const c_PlanarSep* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->getNumberOfPlanes());
}

void c_PlanarSep_getPlane(c_PlanarSep* a_self, const int* a_index,
                          double* a_plane_listed) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0);
  assert(static_cast<IRL::UnsignedIndex_t>(*a_index) <
         a_self->obj_ptr->getNumberOfPlanes());
  a_plane_listed[0] =
      (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
          .normal()[0];
  a_plane_listed[1] =
      (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
          .normal()[1];
  a_plane_listed[2] =
      (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
          .normal()[2];
  a_plane_listed[3] =
      (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
          .distance();
}

bool c_PlanarSep_isFlipped(const c_PlanarSep* a_self) {
  assert(a_self != nullptr);
  return a_self->obj_ptr->isFlipped();
}

void c_PlanarSep_printToScreen(const c_PlanarSep* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  std::cout << (*a_self->obj_ptr);
}

}  // end extern C
