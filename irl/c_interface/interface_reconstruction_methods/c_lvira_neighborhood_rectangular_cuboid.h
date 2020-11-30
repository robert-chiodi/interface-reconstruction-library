// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_LVIRA_NEIGHBORHOOD_RECTANGULAR_CUBOID_H_
#define IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_LVIRA_NEIGHBORHOOD_RECTANGULAR_CUBOID_H_

#include "irl/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "irl/interface_reconstruction_methods/lvira_neighborhood.h"

extern "C" {

struct c_LVIRANeigh_RectCub {
  IRL::LVIRANeighborhood<IRL::RectangularCuboid>* obj_ptr = nullptr;
};

void c_LVIRANeigh_RectCub_new(c_LVIRANeigh_RectCub* a_self);

void c_LVIRANeigh_RectCub_delete(c_LVIRANeigh_RectCub* a_self);

void c_LVIRANeigh_RectCub_setSize(c_LVIRANeigh_RectCub* a_self,
                                  const int* a_size);

void c_LVIRANeigh_RectCub_setMember(c_LVIRANeigh_RectCub* a_self,
                                    const int* a_index,
                                    const c_RectCub* a_rectangular_cuboid,
                                    const double* a_liquid_volume_fraction);

void c_LVIRANeigh_RectCub_addMember(c_LVIRANeigh_RectCub* a_self,
                                    const c_RectCub* a_rectangular_cuboid,
                                    const double* a_volume_fraction);

void c_LVIRANeigh_RectCub_emptyNeighborhood(c_LVIRANeigh_RectCub* a_self);

void c_LVIRANeigh_RectCub_setCenterOfStencil(c_LVIRANeigh_RectCub* a_self,
                                             const int* a_center_cell_index);
}

#endif // IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_LVIRA_NEIGHBORHOOD_RECTANGULAR_CUBOID_H_
