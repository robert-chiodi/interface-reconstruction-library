// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/moments/c_tagged_accumulated_listed_volume_moments_and_normal.h"

extern "C" {

void c_TagAccListVM_VMAN_new(c_TagAccListVM_VMAN* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);

  a_self->obj_ptr = new IRL::TaggedAccumulatedListedVolumeMoments<
      IRL::VolumeMomentsAndNormal>;
}

void c_TagAccListVM_VMAN_delete(c_TagAccListVM_VMAN* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_TagAccListVM_VMAN_getListAtIndex(c_TagAccListVM_VMAN* a_self,
                                        const int* a_list_index,
                                        c_ListVM_VMAN* a_gotten_list) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_gotten_list != nullptr);
  assert(a_gotten_list->obj_ptr != nullptr);
  assert(*a_list_index >= 0);
  assert(*a_list_index < static_cast<int>(a_self->obj_ptr->size()));
  (*a_gotten_list->obj_ptr) =
      (*a_self->obj_ptr)
          .getMomentsForIndex(static_cast<IRL::UnsignedIndex_t>(*a_list_index));
}

void c_TagAccListVM_VMAN_append(c_TagAccListVM_VMAN* a_self,
                                const c_TagAccListVM_VMAN* a_other_list) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_other_list != nullptr);
  (*a_self->obj_ptr) += (*a_other_list->obj_ptr);
}

void c_TagAccListVM_VMAN_clear(c_TagAccListVM_VMAN* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->clear();
}

int c_TagAccListVM_VMAN_getSize(const c_TagAccListVM_VMAN* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->size());
}

int c_TagAccListVM_VMAN_getTagForIndex(const c_TagAccListVM_VMAN* a_self,
                                       const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->getTagForIndex(
      static_cast<IRL::UnsignedIndex_t>(*a_index)));
}
}
