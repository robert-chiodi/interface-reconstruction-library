// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2021 Austin Han <austinhhan@outlook.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/c_interface/interface_reconstruction_methods/c_optimization_behavior.h"

#include <cassert>

extern "C" {

void c_OptimizationBehavior_new(c_OptimizationBehavior* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::OptimizationBehavior;
}

void c_OptimizationBehavior_delete(c_OptimizationBehavior* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

void c_OptimizationBehavior_setAcceptableError(c_OptimizationBehavior* a_self, const double a_parameter) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).acceptable_error = a_parameter;
}

void c_OptimizationBehavior_setMaxIterations(c_OptimizationBehavior* a_self, const int a_parameter) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_parameter >= 0);
  (*a_self->obj_ptr).maximum_iterations = static_cast<IRL::UnsignedIndex_t>(a_parameter);
}

void c_OptimizationBehavior_setMinAngleChange(c_OptimizationBehavior* a_self, const double a_parameter) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).minimum_angle_change = a_parameter;
}

void c_OptimizationBehavior_setMinDistanceChange(c_OptimizationBehavior* a_self, const double a_parameter) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).minimum_distance_change = a_parameter;
}

void c_OptimizationBehavior_setLambdaIncrease(c_OptimizationBehavior* a_self, const double a_parameter) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).lambda_increase = a_parameter;
}

void c_OptimizationBehavior_setLambdaDecrease(c_OptimizationBehavior* a_self, const double a_parameter) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).lambda_decrease = a_parameter;
}

void c_OptimizationBehavior_setDelayJacobianAmount(c_OptimizationBehavior* a_self, const int a_parameter) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_parameter >= 0);
  (*a_self->obj_ptr).delay_jacobian_amount = static_cast<IRL::UnsignedIndex_t>(a_parameter);
}

void c_OptimizationBehavior_setInitialAngle(c_OptimizationBehavior* a_self, const double a_parameter) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).initial_angle = a_parameter;
}

void c_OptimizationBehavior_setInitialDistance(c_OptimizationBehavior* a_self, const double a_parameter) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).initial_distance = a_parameter;
}

void c_OptimizationBehavior_setFiniteDifferenceAngle(c_OptimizationBehavior* a_self, const double a_parameter) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  (*a_self->obj_ptr).finite_difference_angle = a_parameter;
}

void c_OptimizationBehavior_getParameters(c_OptimizationBehavior* a_self, double* a_parameters) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  assert(a_parameters != nullptr);
  a_parameters[0] = (*a_self->obj_ptr).acceptable_error;
  a_parameters[1] = static_cast<double>((*a_self->obj_ptr).maximum_iterations);
  a_parameters[2] = (*a_self->obj_ptr).minimum_angle_change;
  a_parameters[3] = (*a_self->obj_ptr).minimum_distance_change;
  a_parameters[4] = (*a_self->obj_ptr).lambda_increase;
  a_parameters[5] = (*a_self->obj_ptr).lambda_decrease;
  a_parameters[6] = static_cast<double>((*a_self->obj_ptr).delay_jacobian_amount);
  a_parameters[7] = (*a_self->obj_ptr).initial_angle;
  a_parameters[8] = (*a_self->obj_ptr).initial_distance;
  a_parameters[9] = (*a_self->obj_ptr).finite_difference_angle;
}

}  // end extern C
