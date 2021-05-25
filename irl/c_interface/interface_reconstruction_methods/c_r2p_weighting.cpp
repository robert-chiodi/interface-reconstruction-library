// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2021 Austin Han <austinhhan@outlook.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/interface_reconstruction_methods/c_r2p_weighting.h"

extern "C" {

void c_R2PWeighting_new(c_R2PWeighting* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::R2PWeighting;
}

void c_R2PWeighting_delete(c_R2PWeighting* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_R2PWeighting_setImportances(c_R2PWeighting* a_self, const double* a_importances) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).importance_of_liquid_volume_fraction = a_importances[0];
  (*a_self->obj_ptr).importance_of_liquid_centroid_relative_to_gas = a_importances[1];
  (*a_self->obj_ptr).importance_of_centroid = a_importances[2];
  (*a_self->obj_ptr).importance_of_surface_area = a_importances[3];
}

void c_R2PWeighting_setImportanceOfLiquidVolumeFraction(c_R2PWeighting* a_self, const double a_importance) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).importance_of_liquid_volume_fraction = a_importance;
}

void c_R2PWeighting_setImportanceOfLiquidCentroidRelativeToGas(c_R2PWeighting* a_self, const double a_importance) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).importance_of_liquid_centroid_relative_to_gas = a_importance;
}

void c_R2PWeighting_setImportanceOfCentroid(c_R2PWeighting* a_self, const double a_importance) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).importance_of_centroid = a_importance;
}

void c_R2PWeighting_setImportanceOfSurfaceArea(c_R2PWeighting* a_self, const double a_importance) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).importance_of_surface_area = a_importance;
}

void c_R2PWeighting_getImportances(c_R2PWeighting* a_self, double* a_importances) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_importances != nullptr);
  a_importances[0] = (*a_self->obj_ptr).importance_of_liquid_volume_fraction;
  a_importances[1] = (*a_self->obj_ptr).importance_of_liquid_centroid_relative_to_gas;
  a_importances[2] = (*a_self->obj_ptr).importance_of_centroid;
  a_importances[3] = (*a_self->obj_ptr).importance_of_surface_area;
}

}  // end extern C
