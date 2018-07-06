// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_MOMENTS_C_SEPARATED_VOLUME_MOMENTS_H_
#define SRC_C_INTERFACE_MOMENTS_C_SEPARATED_VOLUME_MOMENTS_H_

#include "src/moments/separated_volume_moments.h"

extern "C" {

struct c_SepVM {
  IRL::SeparatedMoments<IRL::VolumeMoments>* obj_ptr = nullptr;
  bool is_owning = false;
};

void c_SepVM_new(c_SepVM* a_self);

void c_SepVM_delete(c_SepVM* a_self);

void c_SepVM_construct(c_SepVM* a_self, const double* a_moments_list);

void c_SepVM_normalizeByVolume(c_SepVM* a_self);

void c_SepVM_multiplyByVolume(c_SepVM* a_self);

double c_SepVM_getVolume(const c_SepVM* a_self, const int* a_index);

void c_SepVM_getCentroid(const c_SepVM* a_self, const int* a_index,
                         double* a_centroid);

const double* c_SepVM_getVolumePtr(const c_SepVM* a_self, const int* a_index);

const double* c_SepVM_getCentroidPtr(const c_SepVM* a_self, const int* a_index);
}

#endif  // SRC_C_INTERFACE_MOMENTS_C_SEPARATED_VOLUME_MOMENTS_H_
