// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2021 Austin Han <austinhhan@outlook.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_OPTIMIZATION_BEHAVIOR_H_
#define IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_OPTIMIZATION_BEHAVIOR_H_

#include "irl/interface_reconstruction_methods/optimization_behavior.h"

extern "C" {

struct c_OptimizationBehavior {
  IRL::OptimizationBehavior* obj_ptr = nullptr;
};

void c_OptimizationBehavior_new(c_OptimizationBehavior* a_self);

void c_OptimizationBehavior_delete(c_OptimizationBehavior* a_self);

void c_OptimizationBehavior_setAcceptableError(c_OptimizationBehavior* a_self, const double a_parameter);

void c_OptimizationBehavior_setMaxIterations(c_OptimizationBehavior* a_self, const int a_parameter);

void c_OptimizationBehavior_setMinAngleChange(c_OptimizationBehavior* a_self, const double a_parameter);

void c_OptimizationBehavior_setMinDistChange(c_OptimizationBehavior* a_self, const double a_parameter);

void c_OptimizationBehavior_setLambdaIncrease(c_OptimizationBehavior* a_self, const double a_parameter);

void c_OptimizationBehavior_setLambdaDecrease(c_OptimizationBehavior* a_self, const double a_parameter);

void c_OptimizationBehavior_setDelayJacobianAmount(c_OptimizationBehavior* a_self, const int a_parameter);

void c_OptimizationBehavior_setInitialAngle(c_OptimizationBehavior* a_self, const double a_parameter);

void c_OptimizationBehavior_setInitialDistance(c_OptimizationBehavior* a_self, const double a_parameter);

void c_OptimizationBehavior_setFiniteDifferenceAngle(c_OptimizationBehavior* a_self, const double a_parameter);

void c_OptimizationBehavior_getParameters(c_OptimizationBehavior* a_self, double* a_parameters);

}  // end extern C

#endif  // IRL_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_OPTIMIZATION_BEHAVIOR_H_
