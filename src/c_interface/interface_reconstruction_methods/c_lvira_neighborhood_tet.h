// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_LVIRA_NEIGHBORHOOD_TET_H_
#define SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_LVIRA_NEIGHBORHOOD_TET_H_

#include "src/c_interface/geometry/polyhedrons/c_tet.h"
#include "src/interface_reconstruction_methods/lvira_neighborhood.h"

extern "C" {

struct c_LVIRANeigh_Tet {
  IRL::LVIRANeighborhood<IRL::Tet>* obj_ptr = nullptr;
};

void c_LVIRANeigh_Tet_new(c_LVIRANeigh_Tet* a_self);

void c_LVIRANeigh_Tet_delete(c_LVIRANeigh_Tet* a_self);

void c_LVIRANeigh_Tet_setSize(c_LVIRANeigh_Tet* a_self,
                                  const int* a_size);

void c_LVIRANeigh_Tet_setMember(c_LVIRANeigh_Tet* a_self,
                                    const int* a_index,
                                    const c_Tet* a_tet,
                                    const double* a_liquid_volume_fraction);

void c_LVIRANeigh_Tet_addMember(c_LVIRANeigh_Tet* a_self,
                                    const c_Tet* a_tet,
                                    const double* a_volume_fraction);

void c_LVIRANeigh_Tet_emptyNeighborhood(c_LVIRANeigh_Tet* a_self);

void c_LVIRANeigh_Tet_setCenterOfStencil(c_LVIRANeigh_Tet* a_self,
                                             const int* a_center_cell_index);
}

#endif  // SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_LVIRA_NEIGHBORHOOD_TET_H_
