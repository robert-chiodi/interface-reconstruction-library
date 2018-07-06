// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_ELVIRA_NEIGHBORHOOD_H_
#define SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_ELVIRA_NEIGHBORHOOD_H_

#include "src/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "src/interface_reconstruction_methods/elvira_neighborhood.h"

extern "C" {

struct c_ELVIRANeigh {
  IRL::ELVIRANeighborhood* obj_ptr = nullptr;
};

void c_ELVIRANeigh_new(c_ELVIRANeigh* a_self);
void c_ELVIRANeigh_delete(c_ELVIRANeigh* a_self);
void c_ELVIRANeigh_setSize(c_ELVIRANeigh* a_self, const int* a_size);
void c_ELVIRANeigh_setMember(c_ELVIRANeigh* a_self,
                             const c_RectCub* a_rectangular_cuboid,
                             const double* a_liquid_volume_fraction,
                             const int* i, const int* j, const int* k);
}

#endif  // SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_ELVIRA_NEIGHBORHOOD_H_
