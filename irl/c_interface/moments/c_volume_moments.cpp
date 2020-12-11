// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/moments/c_volume_moments.h"

extern "C" {

void c_VM_new(c_VM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::VolumeMoments;
  a_self->is_owning = true;
}

void c_VM_delete(c_VM* a_self) {
  if (a_self->is_owning) {
    delete a_self->obj_ptr;
  }
  a_self->obj_ptr = nullptr;
  a_self->is_owning = false;
}

void c_VM_construct(c_VM* a_self, const double* a_moments_list) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_self->is_owning);
  (*a_self->obj_ptr) = IRL::VolumeMoments::fromRawDoublePointer(a_moments_list);
}

void c_VM_normalizeByVolume(c_VM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->normalizeByVolume();
}

void c_VM_multiplyByVolume(c_VM* a_self) {
  assert(a_self != nullptr);
  a_self->obj_ptr->multiplyByVolume();
}

double c_VM_getVolume(const c_VM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return (*a_self->obj_ptr).volume();
}

void c_VM_getCentroid(const c_VM* a_self, double* a_centroid) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_centroid[0] = (*a_self->obj_ptr).centroid()[0];
  a_centroid[1] = (*a_self->obj_ptr).centroid()[1];
  a_centroid[2] = (*a_self->obj_ptr).centroid()[2];
}

const double* c_VM_getVolumePtr(const c_VM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return reinterpret_cast<const double*>(&(*a_self->obj_ptr).volume());
}

const double* c_VM_getCentroidPtr(const c_VM* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return &(*a_self->obj_ptr).centroid()[0];
}
}
