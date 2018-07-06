// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/moments/c_listedvm_vman.h"

extern "C" {

void c_ListVM_VMAN_new(c_ListVM_VMAN* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::ListedVolumeMoments<IRL::VolumeMomentsAndNormal>;
}

void c_ListVM_VMAN_delete(c_ListVM_VMAN* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_ListVM_VMAN_append(c_ListVM_VMAN* a_self, c_ListVM_VMAN* a_other_list) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_other_list != nullptr);
  assert(a_other_list->obj_ptr != nullptr);
  (*a_self->obj_ptr) += (*a_other_list->obj_ptr);
}

void c_ListVM_VMAN_zeroNormalComponent(c_ListVM_VMAN* a_self,
                                       const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);

  for (auto& member : *(a_self->obj_ptr)) {
    member.normalizeByVolume();
    member.normal()[static_cast<IRL::UnsignedIndex_t>(*a_index)] = 0.0;
    member.normal().normalize();
    member.multiplyByVolume();
  }
}

void c_ListVM_VMAN_clear(c_ListVM_VMAN* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->clear();
}

int c_ListVM_VMAN_getSize(const c_ListVM_VMAN* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->size());
}

void c_ListVM_VMAN_getMoments(c_ListVM_VMAN* a_self, const int* a_index,
                              c_VMAN* a_moments) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_moments != nullptr);
  assert(a_moments->obj_ptr != nullptr);
  assert(*a_index >= 0);
  assert(*a_index < static_cast<int>(a_self->obj_ptr->size()));
  (*a_moments->obj_ptr) =
      (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>((*a_index))];
}

void c_ListVM_VMAN_erase(c_ListVM_VMAN* a_self, const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0);
  a_self->obj_ptr->erase(static_cast<IRL::UnsignedIndex_t>((*a_index)));
}
}
