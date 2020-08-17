// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_MOMENTS_C_LISTEDVM_VMAN_H_
#define SRC_C_INTERFACE_MOMENTS_C_LISTEDVM_VMAN_H_

#include "src/c_interface/moments/c_volume_moments_and_normal.h"
#include "src/moments/listed_volume_moments.h"
#include "src/moments/volume_moments_and_normal.h"

extern "C" {

struct c_ListVM_VMAN {
  IRL::ListedVolumeMoments<IRL::VolumeMomentsAndNormal>* obj_ptr = nullptr;
};

void c_ListVM_VMAN_new(c_ListVM_VMAN* a_self);

void c_ListVM_VMAN_delete(c_ListVM_VMAN* a_self);

void c_ListVM_VMAN_append(c_ListVM_VMAN* a_self, c_ListVM_VMAN* a_other_list);

void c_ListVM_VMAN_clear(c_ListVM_VMAN* a_self);

int c_ListVM_VMAN_getSize(const c_ListVM_VMAN* a_self);

void c_ListVM_VMAN_getMoments(c_ListVM_VMAN* a_self, const int* a_index,
                              c_VMAN* a_moments);

void c_ListVM_VMAN_zeroNormalComponent(c_ListVM_VMAN* a_self,
                                       const int* a_index);

void c_ListVM_VMAN_erase(c_ListVM_VMAN* a_self, const int* a_index);
}

#endif  // SRC_C_INTERFACE_MOMENTS_C_LISTEDVM_VMAN_H_
