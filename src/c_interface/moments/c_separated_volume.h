// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_MOMENTS_C_SEPARATED_VOLUME_H_
#define SRC_C_INTERFACE_MOMENTS_C_SEPARATED_VOLUME_H_

#include "src/moments/separated_volume_moments.h"

extern "C" {

struct c_SepVol {
  IRL::SeparatedMoments<IRL::Volume>* obj_ptr = nullptr;
};

void c_SepVol_new(c_SepVol* a_self);

void c_SepVol_delete(c_SepVol* a_self);

void c_SepVol_construct(c_SepVol* a_self, const double* a_moments_list);

double c_SepVol_getVolume(const c_SepVol* a_self, const int* a_index);

}

#endif  // SRC_C_INTERFACE_MOMENTS_C_SEPARATED_VOLUME_H_
