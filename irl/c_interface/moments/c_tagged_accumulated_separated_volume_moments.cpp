// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/moments/c_tagged_accumulated_separated_volume_moments.h"

extern "C" {

void c_TagAccVM_SepVM_new(c_TagAccVM_SepVM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::TaggedAccumulatedVolumeMoments<
      IRL::SeparatedMoments<IRL::VolumeMoments>>;
}

void c_TagAccVM_SepVM_delete(c_TagAccVM_SepVM* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_TagAccVM_SepVM_clear(c_TagAccVM_SepVM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr->clear();
}

void c_TagAccVM_SepVM_normalizeByVolume(c_TagAccVM_SepVM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->normalizeByVolume();
}

void c_TagAccVM_SepVM_multiplyByVolume(c_TagAccVM_SepVM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->multiplyByVolume();
}

void c_TagAccVM_SepVM_getSepVMAtIndex(const c_TagAccVM_SepVM* a_self,
                                    const int* a_index, c_SepVM* a_sep_vm) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index < static_cast<int>(a_self->obj_ptr->size()));
  assert(a_sep_vm != nullptr);
  assert(a_sep_vm->is_owning == false);
  a_sep_vm->obj_ptr = &a_self->obj_ptr->getMomentsForIndex(static_cast<IRL::UnsignedIndex_t>(*a_index));
}

void c_TagAccVM_SepVM_getSepVMAtTag(const c_TagAccVM_SepVM* a_self,
                                    const int* a_tag, c_SepVM* a_sep_vm) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_sep_vm != nullptr);
  assert(a_sep_vm->is_owning == false);
  a_sep_vm->obj_ptr = &a_self->obj_ptr->operator[](static_cast<IRL::UnsignedIndex_t>(*a_tag));
}

int c_TagAccVM_SepVM_getSize(const c_TagAccVM_SepVM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->size());
}

int c_TagAccVM_SepVM_getTagForIndex(const c_TagAccVM_SepVM* a_self,
                                    const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index < static_cast<int>(a_self->obj_ptr->size()));
  return static_cast<int>(a_self->obj_ptr->getTagForIndex(
      static_cast<IRL::UnsignedIndex_t>(*a_index)));
}
}
