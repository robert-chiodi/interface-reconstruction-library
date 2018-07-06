// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/planar_reconstruction/c_planar_separator_path.h"

#include <cassert>
#include <iostream>

extern "C" {

void c_PlanarSepPath_new(c_PlanarSepPath* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::PlanarSeparatorPath();
}

void c_PlanarSepPath_new_PlanarSep(c_PlanarSepPath* a_self, c_PlanarSep* a_planar_separator) {
  assert(a_self->obj_ptr == nullptr);
  assert(a_planar_separator != nullptr);
  assert(a_planar_separator->obj_ptr != nullptr);
  a_self->obj_ptr = new IRL::PlanarSeparatorPath(a_planar_separator->obj_ptr);
}

void c_PlanarSepPath_construct(c_PlanarSepPath* a_self, c_PlanarSep* a_planar_separator) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_planar_separator != nullptr);
  assert(a_planar_separator->obj_ptr != nullptr);
  *(a_self->obj_ptr) = IRL::PlanarSeparatorPath(a_planar_separator->obj_ptr);
}

void c_PlanarSepPath_delete(c_PlanarSepPath* a_self) {
  a_self->obj_ptr = nullptr;
}

void c_PlanarSepPath_setEdgeConnectivity(c_PlanarSepPath* a_self,
                                   const c_PlanarSepPath* a_self_to_neighbor) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_self_to_neighbor != nullptr);
  a_self->obj_ptr->setEdgeConnectivity(
      a_self_to_neighbor->obj_ptr);
}

void c_PlanarSepPath_setEdgeConnectivityNull(c_PlanarSepPath* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->setEdgeConnectivity(nullptr);
}

void c_PlanarSepPath_setId(c_PlanarSepPath* a_self, const int* a_id) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->setId(static_cast<IRL::UnsignedIndex_t>(*a_id));
}

int c_PlanarSepPath_getId(c_PlanarSepPath* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->getId());
}

} // end extern "C"
