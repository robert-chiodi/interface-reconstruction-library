// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/moments/c_separated_volume_moments.h"

extern "C" {

void c_SepVM_new(c_SepVM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::SeparatedMoments<IRL::VolumeMoments>;
  a_self->is_owning = true;
}

void c_SepVM_delete(c_SepVM* a_self) {
  if(a_self->is_owning){
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}

void c_SepVM_construct(c_SepVM* a_self, const double* a_moments_list) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_self->is_owning);
  (*a_self->obj_ptr) =
      IRL::SeparatedMoments<IRL::VolumeMoments>::fromRawDoublePointer(
          a_moments_list);
}

void c_SepVM_normalizeByVolume(c_SepVM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->normalizeByVolume();
}

void c_SepVM_multiplyByVolume(c_SepVM* a_self) {
  assert(a_self != nullptr);
  a_self->obj_ptr->multiplyByVolume();
}

double c_SepVM_getVolume(const c_SepVM* a_self, const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0);
  assert(*a_index < 2);
  return (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
      .volume();
}

void c_SepVM_getCentroid(const c_SepVM* a_self, const int* a_index,
                         double* a_centroid) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0);
  assert(*a_index < 2);
  a_centroid[0] =
      (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
          .centroid()[0];
  a_centroid[1] =
      (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
          .centroid()[1];
  a_centroid[2] =
      (*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
          .centroid()[2];
}

const double* c_SepVM_getVolumePtr(const c_SepVM* a_self, const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0);
  assert(*a_index < 2);
  return reinterpret_cast<const double*>(
      &(*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
           .volume());
}

const double* c_SepVM_getCentroidPtr(const c_SepVM* a_self,
                                     const int* a_index) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(*a_index >= 0);
  assert(*a_index < 2);
  return &(*a_self->obj_ptr)[static_cast<IRL::UnsignedIndex_t>(*a_index)]
              .centroid()[0];
}
}
