// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_TAGGED_ACCUMULATED_VOLUME_H_
#define SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_TAGGED_ACCUMULATED_VOLUME_H_

#include "src/moments/volume.h"
#include "src/moments/tagged_accumulated_volume_moments.h"
#include "src/c_interface/moments/c_tagged_accumulated_volume.h"

extern "C" {

struct c_TagAccVM2_Vol {
  IRL::TaggedAccumulatedVolumeMoments<IRL::TaggedAccumulatedVolumeMoments<
      IRL::Volume>>* obj_ptr = nullptr;
  bool is_owning = false;
};

void c_TagAccVM2_Vol_new(c_TagAccVM2_Vol* a_self);

void c_TagAccVM2_Vol_delete(c_TagAccVM2_Vol* a_self);

void c_TagAccVM2_Vol_getAtIndex(const c_TagAccVM2_Vol* a_self,
                                const int* a_list_index, c_TagAccVM_Vol* a_tagged_volume);

void c_TagAccVM2_Vol_getAtTag(const c_TagAccVM2_Vol* a_self,
                              const int* a_tag, c_TagAccVM_Vol* a_tagged_volume);

int c_TagAccVM2_Vol_getSize(const c_TagAccVM2_Vol* a_self);

int c_TagAccVM2_Vol_getTagForIndex(const c_TagAccVM2_Vol* a_self,
                                    const int* a_index);
}

#endif  // SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_TAGGED_ACCUMULATED_VOLUME_H_
