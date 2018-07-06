// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/moments/c_volume_moments_and_normal.h"

#include <cassert>

extern "C" {

void c_VMAN_new(c_VMAN* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::VolumeMomentsAndNormal;
}

void c_VMAN_delete(c_VMAN* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_VMAN_getVolume(c_VMAN* a_self, double* a_volume) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  *a_volume = a_self->obj_ptr->volumeMoments().volume();
}

void c_VMAN_getCentroid(c_VMAN* a_self, double* a_centroid) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_centroid[0] = a_self->obj_ptr->volumeMoments().centroid()[0];
  a_centroid[1] = a_self->obj_ptr->volumeMoments().centroid()[1];
  a_centroid[2] = a_self->obj_ptr->volumeMoments().centroid()[2];
}

void c_VMAN_getNormal(c_VMAN* a_self, double* a_normal) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_normal[0] = a_self->obj_ptr->normal()[0];
  a_normal[1] = a_self->obj_ptr->normal()[1];
  a_normal[2] = a_self->obj_ptr->normal()[2];
}

void c_VMAN_normalizeByVolume(c_VMAN* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->normalizeByVolume();
}

void c_VMAN_multiplyByVolume(c_VMAN* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->multiplyByVolume();
}
}
