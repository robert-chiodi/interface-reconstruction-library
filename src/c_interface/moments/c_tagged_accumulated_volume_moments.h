// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_VOLUME_MOMENTS_H_
#define SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_VOLUME_MOMENTS_H_

#include "src/moments/tagged_accumulated_volume_moments.h"
#include "src/moments/volume_moments.h"

extern "C" {

struct c_TagAccVM_VM {
  IRL::TaggedAccumulatedVolumeMoments<IRL::VolumeMoments>* obj_ptr = nullptr;
};

void c_TagAccVM_VM_new(c_TagAccVM_VM* a_self);

void c_TagAccVM_VM_delete(c_TagAccVM_VM* a_self);

void c_TagAccVM_VM_normalizeByVolume(c_TagAccVM_VM* a_self);

void c_TagAccVM_VM_multiplyByVolume(c_TagAccVM_VM* a_self);

double c_TagAccVM_VM_getVolumeAtIndex(const c_TagAccVM_VM* a_self,
                                      const int* a_list_index);

void c_TagAccVM_VM_getCentroidAtIndex(c_TagAccVM_VM* a_self,
                                      const int* a_list_index,
                                      double* a_centroid);

const double* c_TagAccVM_VM_getVolumePtrAtIndex(const c_TagAccVM_VM* a_self,
                                                const int* a_list_index);

const double* c_TagAccVM_VM_getCentroidPtrAtIndex(c_TagAccVM_VM* a_self,
                                                  const int* a_list_index);

int c_TagAccVM_VM_getSize(const c_TagAccVM_VM* a_self);

int c_TagAccVM_VM_getTagForIndex(const c_TagAccVM_VM* a_self,
                                 const int* a_index);
}

#endif  // SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_VOLUME_MOMENTS_H_
