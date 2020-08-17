// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_MOMENTS_C_SEPARATED_VOLUME_MOMENTS_DOUBLES3_H_
#define SRC_C_INTERFACE_MOMENTS_C_SEPARATED_VOLUME_MOMENTS_DOUBLES3_H_

#include "src/moments/separated_volume_moments.h"
#include "src/moments/volume_moments_and_doubles.h"

extern "C" {

struct c_SepVM_d3 {
  IRL::SeparatedMoments<IRL::VolumeMomentsAndDoubles<3>>* obj_ptr = nullptr;
};

void c_SepVM_d3_new(c_SepVM_d3* a_self);

void c_SepVM_d3_delete(c_SepVM_d3* a_self);

void c_SepVM_d3_normalizeByVolume(c_SepVM_d3* a_self);

void c_SepVM_d3_multiplyByVolume(c_SepVM_d3* a_self);

double c_SepVM_d3_getVolume(const c_SepVM_d3* a_self, const int* a_index);

void c_SepVM_d3_getCentroid(const c_SepVM_d3* a_self, const int* a_index,
                            double* a_centroid);

void c_SepVM_d3_getData(const c_SepVM_d3* a_self, const int* a_index,
                        double* a_data);

const double* c_SepVM_d3_getVolumePtr(const c_SepVM_d3* a_self,
                                      const int* a_index);

const double* c_SepVM_d3_getCentroidPtr(const c_SepVM_d3* a_self,
                                        const int* a_index);
}

#endif  // SRC_C_INTERFACE_MOMENTS_C_SEPARATED_VOLUME_MOMENTS_DOUBLES3_H_
