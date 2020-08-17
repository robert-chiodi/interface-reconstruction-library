// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/interface_reconstruction_methods/c_elvira_neighborhood.h"

#include <cassert>

extern "C" {

void c_ELVIRANeigh_new(c_ELVIRANeigh* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::ELVIRANeighborhood;
}

void c_ELVIRANeigh_delete(c_ELVIRANeigh* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_ELVIRANeigh_setSize(c_ELVIRANeigh* a_self, const int* a_size) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->resize(static_cast<IRL::UnsignedIndex_t>(*a_size));
}

void c_ELVIRANeigh_setMember(c_ELVIRANeigh* a_self,
                             const c_RectCub* a_rectangular_cuboid,
                             const double* a_liquid_volume_fraction,
                             const int* i, const int* j, const int* k) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_rectangular_cuboid != nullptr);
  assert(a_rectangular_cuboid->obj_ptr != nullptr);
  a_self->obj_ptr->setMember(a_rectangular_cuboid->obj_ptr,
                             a_liquid_volume_fraction, *i, *j, *k);
}
}
