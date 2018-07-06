// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/interface_reconstruction_methods/c_lvira_neighborhood_hexahedron.h"

#include <cassert>

extern "C" {

void c_LVIRANeigh_Hex_new(c_LVIRANeigh_Hex* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::LVIRANeighborhood<IRL::Hexahedron>;
}

void c_LVIRANeigh_Hex_delete(c_LVIRANeigh_Hex* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_LVIRANeigh_Hex_setSize(c_LVIRANeigh_Hex* a_self,
                                  const int* a_size) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->resize(static_cast<IRL::UnsignedIndex_t>(*a_size));
}

void c_LVIRANeigh_Hex_setMember(c_LVIRANeigh_Hex* a_self,
                                    const int* a_index,
                                    const c_Hex* a_hexahedron,
                                    const double* a_liquid_volume_fraction) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_hexahedron != nullptr);
  assert(a_hexahedron->obj_ptr != nullptr);
  assert(*a_index >= 0);
  assert(*a_index < static_cast<int>(a_self->obj_ptr->size()));
  a_self->obj_ptr->setMember(static_cast<IRL::UnsignedIndex_t>(*a_index),
		  a_hexahedron->obj_ptr,
                             a_liquid_volume_fraction);
}

void c_LVIRANeigh_Hex_addMember(c_LVIRANeigh_Hex* a_self,
                                    const c_Hex* a_hexahedron,
                                    const double* a_volume_fraction) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_hexahedron != nullptr);
  assert(a_hexahedron->obj_ptr != nullptr);
  assert(a_volume_fraction != nullptr);
  a_self->obj_ptr->addMember(a_hexahedron->obj_ptr, a_volume_fraction);
}

void c_LVIRANeigh_Hex_emptyNeighborhood(c_LVIRANeigh_Hex* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->emptyNeighborhood();
}

void c_LVIRANeigh_Hex_setCenterOfStencil(c_LVIRANeigh_Hex* a_self,
                                             const int* a_center_cell_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_center_cell_index >= 0);
  assert(*a_center_cell_index < static_cast<int>(a_self->obj_ptr->size()));
  a_self->obj_ptr->setCenterOfStencil(
      static_cast<IRL::UnsignedIndex_t>(*a_center_cell_index));
}

}  // end extern "C"
