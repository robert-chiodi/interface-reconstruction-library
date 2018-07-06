// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/interface_reconstruction_methods/c_r2p_neighborhood_rectangular_cuboid.h"

extern "C" {

void c_R2PNeigh_RectCub_new(c_R2PNeigh_RectCub* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::R2PNeighborhood<IRL::RectangularCuboid>;
}

void c_R2PNeigh_RectCub_delete(c_R2PNeigh_RectCub* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_R2PNeigh_RectCub_setSize(c_R2PNeigh_RectCub* a_self, const int* a_size) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_size >= 0);
  a_self->obj_ptr->resize(static_cast<IRL::UnsignedIndex_t>(*a_size));
}

void c_R2PNeigh_RectCub_setMember(c_R2PNeigh_RectCub* a_self,
                                  const c_RectCub* a_rectangular_cuboid,
                                  const c_SepVM* a_separated_volume_moments,
                                  const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_separated_volume_moments != nullptr);
  assert(a_separated_volume_moments->obj_ptr != nullptr);
  assert(*a_index >= 0);
  a_self->obj_ptr->setMember(static_cast<IRL::UnsignedIndex_t>(*a_index),
                             a_rectangular_cuboid->obj_ptr,
                             a_separated_volume_moments->obj_ptr);
}

void c_R2PNeigh_RectCub_addMember(c_R2PNeigh_RectCub* a_self,
                                  const c_RectCub* a_rectangular_cuboid,
                                  const c_SepVM* a_separated_volume_moments) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  assert(a_separated_volume_moments != nullptr);
  assert(a_separated_volume_moments->obj_ptr != nullptr);
  a_self->obj_ptr->addMember(a_rectangular_cuboid->obj_ptr,
                             a_separated_volume_moments->obj_ptr);
}

void c_R2PNeigh_RectCub_emptyNeighborhood(c_R2PNeigh_RectCub* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->emptyNeighborhood();
}

void c_R2PNeigh_RectCub_setCenterOfStencil(c_R2PNeigh_RectCub* a_self,
                                           const int* a_center_cell_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_center_cell_index >= 0);
  assert(*a_center_cell_index < static_cast<int>(a_self->obj_ptr->size()));
  a_self->obj_ptr->setCenterOfStencil(
      static_cast<IRL::UnsignedIndex_t>(*a_center_cell_index));
}

void c_R2PNeigh_RectCub_setSurfaceArea(c_R2PNeigh_RectCub* a_self,
                                       const double* a_surface_area) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->setSurfaceArea(*a_surface_area);
}

}  // end extern C
