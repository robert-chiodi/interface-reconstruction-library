// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/moments/c_tagged_accumulated_separated_volume.h"

extern "C" {

void c_TagAccVM_SepVol_new(c_TagAccVM_SepVol* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::TaggedAccumulatedVolumeMoments<
      IRL::SeparatedMoments<IRL::Volume>>;
}

void c_TagAccVM_SepVol_delete(c_TagAccVM_SepVol* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}


double c_TagAccVM_SepVol_getVolumeAtIndex(const c_TagAccVM_SepVol* a_self,
                                         const int* a_list_index,
                                         const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_list_index >= 0);
  assert(*a_list_index < static_cast<int>(a_self->obj_ptr->size()));
  assert(*a_index >= 0);
  assert(*a_index < 2);
  return a_self->obj_ptr
      ->getMomentsForIndex(static_cast<IRL::UnsignedIndex_t>(
          *a_list_index))[static_cast<IRL::UnsignedIndex_t>(*a_index)];
}

double c_TagAccVM_SepVol_getVolumeAtTag(const c_TagAccVM_SepVol* a_self,
                                       const int* a_tag, const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_tag >= 0);
  assert(
      a_self->obj_ptr->isTagKnown(static_cast<IRL::UnsignedIndex_t>(*a_tag)));
  assert(*a_index >= 0);
  assert(*a_index < 2);
  return (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_tag)]
                           [static_cast<IRL::UnsignedIndex_t>(*a_index)];
}

const double* c_TagAccVM_SepVol_getVolumePtrAtIndex(
    const c_TagAccVM_SepVol* a_self, const int* a_list_index,
    const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_list_index >= 0);
  assert(*a_list_index < static_cast<int>(a_self->obj_ptr->size()));
  assert(*a_index >= 0);
  assert(*a_index < 2);
  return reinterpret_cast<const double*>(
      &a_self->obj_ptr
           ->getMomentsForIndex(static_cast<IRL::UnsignedIndex_t>(
               *a_list_index))[static_cast<IRL::UnsignedIndex_t>(*a_index)]);
}

int c_TagAccVM_SepVol_getSize(const c_TagAccVM_SepVol* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->size());
}

int c_TagAccVM_SepVol_getTagForIndex(const c_TagAccVM_SepVol* a_self,
                                    const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index < static_cast<int>(a_self->obj_ptr->size()));
  return static_cast<int>(a_self->obj_ptr->getTagForIndex(
      static_cast<IRL::UnsignedIndex_t>(*a_index)));
}
}
