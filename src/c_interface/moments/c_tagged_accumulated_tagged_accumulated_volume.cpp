// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/moments/c_tagged_accumulated_tagged_accumulated_volume.h"

extern "C" {

void c_TagAccVM2_Vol_new(c_TagAccVM2_Vol* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::TaggedAccumulatedVolumeMoments<IRL::TaggedAccumulatedVolumeMoments<
      IRL::Volume>>;
  a_self->is_owning = true;
}

void c_TagAccVM2_Vol_delete(c_TagAccVM2_Vol* a_self) {
  assert(a_self != nullptr);
  if(a_self->is_owning){
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}


void c_TagAccVM2_Vol_getAtIndex(const c_TagAccVM2_Vol* a_self,
                                         const int* a_list_index, c_TagAccVM_Vol* a_tagged_volume){
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_list_index >= 0);
  assert(*a_list_index < static_cast<int>(a_self->obj_ptr->size()));
  assert(a_tagged_volume->is_owning == false);
  a_tagged_volume->obj_ptr = &(a_self->obj_ptr
      ->getMomentsForIndex(static_cast<IRL::UnsignedIndex_t>(
          *a_list_index)));
}

void c_TagAccVM2_Vol_getAtTag(const c_TagAccVM2_Vol* a_self,
                                      const int* a_tag, c_TagAccVM_Vol* a_tagged_volume){
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_tag >= 0);
  assert(
      a_self->obj_ptr->isTagKnown(static_cast<IRL::UnsignedIndex_t>(*a_tag)));
  assert(a_tagged_volume->is_owning == false);
  a_tagged_volume->obj_ptr = &((*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_tag)]);
}

int c_TagAccVM2_Vol_getSize(const c_TagAccVM2_Vol* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->size());
}

int c_TagAccVM2_Vol_getTagForIndex(const c_TagAccVM2_Vol* a_self,
                                    const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index < static_cast<int>(a_self->obj_ptr->size()));
  return static_cast<int>(a_self->obj_ptr->getTagForIndex(
      static_cast<IRL::UnsignedIndex_t>(*a_index)));
}
}
