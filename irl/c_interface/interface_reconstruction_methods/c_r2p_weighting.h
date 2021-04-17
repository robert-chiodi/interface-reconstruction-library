// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2021 Austin Han <austinhhan@outlook.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_R2P_WEIGHTING_H_
#define IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_R2P_WEIGHTING_H_

#include "irl/interface_reconstruction_methods/r2p_optimization.h"

extern "C" {

struct c_R2PWeighting {
  IRL::R2PWeighting* obj_ptr = nullptr;
};

void c_R2PWeighting_new(c_R2PWeighting* a_self);

void c_R2PWeighting_delete(c_R2PWeighting* a_self);

void c_R2PWeighting_setImportances(c_R2PWeighting* a_self, const double* a_importances);

void c_R2PWeighting_setImpOfLiqVolFrac(c_R2PWeighting* a_self, const double* a_importance);

void c_R2PWeighting_setImpOfLiqCentRelToGas(c_R2PWeighting* a_self, const double* a_importance);

void c_R2PWeighting_setImpOfCentroid(c_R2PWeighting* a_self, const double* a_importance);

void c_R2PWeighting_setImpOfSurfArea(c_R2PWeighting* a_self, const double* a_importance);

void c_R2PWeighting_getImportances(c_R2PWeighting* a_self, double* a_importances);

}  // end extern C

#endif  // IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_R2P_WEIGHTING_H_
