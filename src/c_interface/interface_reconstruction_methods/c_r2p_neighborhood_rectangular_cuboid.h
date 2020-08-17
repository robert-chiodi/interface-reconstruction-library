// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_R2P_NEIGHBORHOOD_RECTANGULAR_CUBOID_H_
#define SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_R2P_NEIGHBORHOOD_RECTANGULAR_CUBOID_H_

#include "src/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "src/c_interface/moments/c_separated_volume_moments.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/interface_reconstruction_methods/r2p_optimization.h"

extern "C" {

struct c_R2PNeigh_RectCub {
  IRL::R2PNeighborhood<IRL::RectangularCuboid>* obj_ptr = nullptr;
};

void c_R2PNeigh_RectCub_new(c_R2PNeigh_RectCub* a_self);

void c_R2PNeigh_RectCub_delete(c_R2PNeigh_RectCub* a_self);

void c_R2PNeigh_RectCub_setSize(c_R2PNeigh_RectCub* a_self, const int* a_size);

void c_R2PNeigh_RectCub_setMember(c_R2PNeigh_RectCub* a_self,
                                  const c_RectCub* a_rectangular_cuboid,
                                  const c_SepVM* a_separated_volume_moments,
                                  const int* a_index);

void c_R2PNeigh_RectCub_addMember(c_R2PNeigh_RectCub* a_self,
                                  const c_RectCub* a_rectangular_cuboid,
                                  const c_SepVM* a_separated_volume_moments);

void c_R2PNeigh_RectCub_emptyNeighborhood(c_R2PNeigh_RectCub* a_self);

void c_R2PNeigh_RectCub_setCenterOfStencil(c_R2PNeigh_RectCub* a_self,
                                           const int* a_center_cell_index);

void c_R2PNeigh_RectCub_setSurfaceArea(c_R2PNeigh_RectCub* a_self,
                                       const double* a_surface_area);

}  // end extern C

#endif  // SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_R2P_NEIGHBORHOOD_RECTANGULAR_CUBOID_H_
