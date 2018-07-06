// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_LISTED_VOLUME_MOMENTS_AND_NORMAL_H_
#define SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_LISTED_VOLUME_MOMENTS_AND_NORMAL_H_

#include "src/c_interface/moments/c_listedvm_vman.h"
#include "src/moments/tagged_accumulated_listed_volume_moments.h"
#include "src/moments/volume_moments_and_normal.h"

extern "C" {

struct c_TagAccListVM_VMAN {
  IRL::TaggedAccumulatedListedVolumeMoments<IRL::VolumeMomentsAndNormal>*
      obj_ptr = nullptr;
};

void c_TagAccListVM_VMAN_new(c_TagAccListVM_VMAN* a_self);

void c_TagAccListVM_VMAN_delete(c_TagAccListVM_VMAN* a_self);

void c_TagAccListVM_VMAN_getListAtIndex(c_TagAccListVM_VMAN* a_self,
                                        const int* a_list_index,
                                        c_ListVM_VMAN* a_gotten_list);

void c_TagAccListVM_VMAN_append(c_TagAccListVM_VMAN* a_self,
                                const c_TagAccListVM_VMAN* a_other_list);

void c_TagAccListVM_VMAN_clear(c_TagAccListVM_VMAN* a_self);

int c_TagAccListVM_VMAN_getSize(const c_TagAccListVM_VMAN* a_self);

int c_TagAccListVM_VMAN_getTagForIndex(const c_TagAccListVM_VMAN* a_self,
                                       const int* a_index);
}

#endif  // SRC_C_INTERFACE_MOMENTS_C_TAGGED_ACCUMULATED_LISTED_VOLUME_MOMENTS_AND_NORMAL_H_
