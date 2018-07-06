// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_PLANAR_SEPARATOR_PATH_GROUP_H_
#define SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_PLANAR_SEPARATOR_PATH_GROUP_H_

#include "src/planar_reconstruction/planar_separator_path_group.h"
#include "src/c_interface/planar_reconstruction/c_planar_separator_path.h"
#include "src/parameters/defined_types.h"

extern "C" {

struct c_PlanarSepPathGroup {
  IRL::PlanarSeparatorPathGroup* obj_ptr = nullptr;
};

void c_PlanarSepPathGroup_new(c_PlanarSepPathGroup* a_self);

void c_PlanarSepPathGroup_delete(c_PlanarSepPathGroup* a_self);

void c_PlanarSepPathGroup_addPlanarSeparatorPath(c_PlanarSepPathGroup* a_self,
												 const c_PlanarSepPath* a_planar_separator_path);

void c_PlanarSepPathGroup_addPlanarSeparatorPath_Id(c_PlanarSepPathGroup* a_self,
												 const c_PlanarSepPath* a_planar_separator_path,
												 const int* a_id);

void c_PlanarSepPathGroup_setPriorityOrder(c_PlanarSepPathGroup* a_self,
										   const int* a_size,
										   const int* a_priority_order);

int c_PlanarSepPathGroup_getPriorityOrderSize(c_PlanarSepPathGroup* a_self);

int c_PlanarSepPathGroup_getPriorityOrderTag(c_PlanarSepPathGroup* a_self, const int* a_index);

}
#endif  // SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_PLANAR_SEPARATOR_PATH_GROUP_H_
