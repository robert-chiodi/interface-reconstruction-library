// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_MOMENTS_C_VOLUME_MOMENTS_H_
#define IRL_C_INTERFACE_MOMENTS_C_VOLUME_MOMENTS_H_

#include "irl/moments/volume_moments.h"

extern "C" {

struct c_VM {
  IRL::VolumeMoments* obj_ptr = nullptr;
  bool is_owning = false;
};

void c_VM_new(c_VM* a_self);

void c_VM_delete(c_VM* a_self);

void c_VM_construct(c_VM* a_self, const double* a_moments_list);

void c_VM_normalizeByVolume(c_VM* a_self);

void c_VM_multiplyByVolume(c_VM* a_self);

double c_VM_getVolume(const c_VM* a_self);

void c_VM_getCentroid(const c_VM* a_self, double* a_centroid);

const double* c_VM_getVolumePtr(const c_VM* a_self);

const double* c_VM_getCentroidPtr(const c_VM* a_self);
}

#endif  // IRL_C_INTERFACE_MOMENTS_C_VOLUME_MOMENTS_H_
