// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/moments/c_tagged_accumulated_volume_moments.h"

#include <cassert>
#include <cstring>

extern "C" {

void c_TagAccVM_VM_new(c_TagAccVM_VM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::TaggedAccumulatedVolumeMoments<IRL::VolumeMoments>;
}

void c_TagAccVM_VM_delete(c_TagAccVM_VM* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_TagAccVM_VM_normalizeByVolume(c_TagAccVM_VM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->normalizeByVolume();
}

void c_TagAccVM_VM_multiplyByVolume(c_TagAccVM_VM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->multiplyByVolume();
}

double c_TagAccVM_VM_getVolumeAtIndex(const c_TagAccVM_VM* a_self,
                                      const int* a_list_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_list_index >= 0);
  assert(*a_list_index < static_cast<int>(a_self->obj_ptr->size()));
  return a_self->obj_ptr
      ->getMomentsForIndex(static_cast<IRL::UnsignedIndex_t>(*a_list_index))
      .volume();
}

void c_TagAccVM_VM_getCentroidAtIndex(c_TagAccVM_VM* a_self,
                                      const int* a_list_index,
                                      double* a_centroid) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_list_index >= 0);
  assert(*a_list_index < static_cast<int>(a_self->obj_ptr->size()));
  IRL::Pt centroid =
      a_self->obj_ptr
          ->getMomentsForIndex(static_cast<IRL::UnsignedIndex_t>(*a_list_index))
          .centroid();
  std::memcpy(a_centroid, &centroid[0], 3 * sizeof(double));
}

const double* c_TagAccVM_VM_getVolumePtrAtIndex(const c_TagAccVM_VM* a_self,
                                                const int* a_list_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_list_index >= 0);
  assert(*a_list_index < static_cast<int>(a_self->obj_ptr->size()));
  return reinterpret_cast<const double*>(
      &a_self->obj_ptr
           ->getMomentsForIndex(
               static_cast<IRL::UnsignedIndex_t>(*a_list_index))
           .volume());
}

const double* c_TagAccVM_VM_getCentroidPtrAtIndex(c_TagAccVM_VM* a_self,
                                                  const int* a_list_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_list_index >= 0);
  assert(*a_list_index < static_cast<int>(a_self->obj_ptr->size()));
  return &a_self->obj_ptr
              ->getMomentsForIndex(
                  static_cast<IRL::UnsignedIndex_t>(*a_list_index))
              .centroid()[0];
}

int c_TagAccVM_VM_getSize(const c_TagAccVM_VM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->size());
}

int c_TagAccVM_VM_getTagForIndex(const c_TagAccVM_VM* a_self,
                                 const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return static_cast<int>(a_self->obj_ptr->getTagForIndex(
      static_cast<IRL::UnsignedIndex_t>(*a_index)));
}
}
