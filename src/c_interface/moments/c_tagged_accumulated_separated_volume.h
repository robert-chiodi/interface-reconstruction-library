// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_SEPARATED_VOLUME_H_
#define SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_SEPARATED_VOLUME_H_

#include "src/moments/separated_volume_moments.h"
#include "src/moments/tagged_accumulated_volume_moments.h"

extern "C" {

struct c_TagAccVM_SepVol {
  IRL::TaggedAccumulatedVolumeMoments<
      IRL::SeparatedMoments<IRL::Volume>>* obj_ptr = nullptr;
};

void c_TagAccVM_SepVol_new(c_TagAccVM_SepVol* a_self);

void c_TagAccVM_SepVol_delete(c_TagAccVM_SepVol* a_self);

double c_TagAccVM_SepVol_getVolumeAtIndex(const c_TagAccVM_SepVol* a_self,
                                         const int* a_list_index,
                                         const int* a_index);

double c_TagAccVM_SepVol_getVolumeAtTag(const c_TagAccVM_SepVol* a_self,
                                       const int* a_tag, const int* a_index);

const double* c_TagAccVM_SepVol_getVolumePtrAtIndex(
    const c_TagAccVM_SepVol* a_self, const int* a_list_index,
    const int* a_index);

int c_TagAccVM_SepVol_getSize(const c_TagAccVM_SepVol* a_self);

int c_TagAccVM_SepVol_getTagForIndex(const c_TagAccVM_SepVol* a_self,
                                    const int* a_index);
}

#endif  // SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_SEPARATED_VOLUME_H_
