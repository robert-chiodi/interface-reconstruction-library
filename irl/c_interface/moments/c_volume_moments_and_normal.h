// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_MOMENTS_C_VOLUME_MOMENTS_AND_NORMAL_H_
#define IRL_C_INTERFACE_MOMENTS_C_VOLUME_MOMENTS_AND_NORMAL_H_

#include "irl/moments/volume_moments_and_normal.h"

extern "C" {

struct c_VMAN {
  IRL::VolumeMomentsAndNormal* obj_ptr = nullptr;
};

void c_VMAN_new(c_VMAN* a_self);

void c_VMAN_delete(c_VMAN* a_self);

void c_VMAN_getVolume(c_VMAN* a_self, double* a_volume);

void c_VMAN_getCentroid(c_VMAN* a_self, double* a_centroid);

void c_VMAN_getNormal(c_VMAN* a_self, double* a_normal);

void c_VMAN_normalizeByVolume(c_VMAN* a_self);

void c_VMAN_multiplyByVolume(c_VMAN* a_self);
}

#endif // IRL_C_INTERFACE_MOMENTS_C_VOLUME_MOMENTS_AND_NORMAL_H_
