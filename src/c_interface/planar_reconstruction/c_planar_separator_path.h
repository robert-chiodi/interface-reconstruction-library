// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_PLANAR_SEPARATOR_PATH_H_
#define SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_PLANAR_SEPARATOR_PATH_H_

#include "src/planar_reconstruction/planar_separator_path.h"
#include "src/c_interface/planar_reconstruction/c_separators.h"
#include "src/parameters/defined_types.h"

extern "C" {

struct c_PlanarSepPath {
  IRL::PlanarSeparatorPath* obj_ptr = nullptr;
};

void c_PlanarSepPath_new(c_PlanarSepPath* a_self);

void c_PlanarSepPath_new_PlanarSep(c_PlanarSepPath* a_self, c_PlanarSep* a_planar_separator);

void c_PlanarSepPath_construct(c_PlanarSepPath* a_self, c_PlanarSep* a_planar_separator);

void c_PlanarSepPath_delete(c_PlanarSepPath* a_self);

void c_PlanarSepPath_setEdgeConnectivity(c_PlanarSepPath* a_self,
										 const c_PlanarSepPath* a_ptr_to_neighbor);

void c_PlanarSepPath_setEdgeConnectivityNull(c_PlanarSepPath* a_self);

void c_PlanarSepPath_setId(c_PlanarSepPath* a_self, const int* a_id);

int c_PlanarSepPath_getId(c_PlanarSepPath* a_self);

}
#endif  // SRC_C_INTERFACE_PLANAR_RECONSTRUCTION_C_PLANAR_SEPARATOR_PATH_H_
