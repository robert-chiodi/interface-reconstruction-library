// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_SEPARATED_VOLUME_MOMENTS_H_
#define SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_SEPARATED_VOLUME_MOMENTS_H_

#include "src/moments/separated_volume_moments.h"
#include "src/moments/tagged_accumulated_listed_volume_moments.h"

#include "src/c_interface/moments/c_separated_volume_moments.h"

extern "C" {

struct c_TagAccVM_SepVM {
  IRL::TaggedAccumulatedVolumeMoments<
      IRL::SeparatedMoments<IRL::VolumeMoments>>* obj_ptr = nullptr;
};

void c_TagAccVM_SepVM_new(c_TagAccVM_SepVM* a_self);

void c_TagAccVM_SepVM_delete(c_TagAccVM_SepVM* a_self);

void c_TagAccVM_SepVM_normalizeByVolume(c_TagAccVM_SepVM* a_self);

void c_TagAccVM_SepVM_multiplyByVolume(c_TagAccVM_SepVM* a_self);

void c_TagAccVM_SepVM_getSepVMAtIndex(const c_TagAccVM_SepVM* a_self,
                                      const int* a_index, c_SepVM* a_sep_vm);

void c_TagAccVM_SepVM_getSepVMAtTag(const c_TagAccVM_SepVM* a_self,
                                    const int* a_tag, c_SepVM* a_sep_vm);

int c_TagAccVM_SepVM_getSize(const c_TagAccVM_SepVM* a_self);

int c_TagAccVM_SepVM_getTagForIndex(const c_TagAccVM_SepVM* a_self,
                                    const int* a_index);
}

#endif  // SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_SEPARATED_VOLUME_MOMENTS_H_
