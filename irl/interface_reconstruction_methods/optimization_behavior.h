// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2021 Austin Han <han.austin@outlook.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_OPTIMIZATION_BEHAVIOR_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_OPTIMIZATION_BEHAVIOR_H_

#include "irl/parameters/defined_types.h"

namespace IRL {

/// \brief Struct to contain LVIRA/MoF/R2P optimization parameters
struct OptimizationBehavior {
  /// \brief If `this->calculateScalarError()` is less than this, exit.
  double acceptable_error = 1.0e-4 * 1.0e-4;
  /// \brief Maximum number of attempted iterations before exiting.
  UnsignedIndex_t maximum_iterations = 20;
  /// \brief Minimum change in angle related delta below which minimum is
  /// deemed reached.
  double minimum_angle_change = 0.0001745329;
  /// \brief Minimum change in distance related delta below which minimum is
  /// deemed reached.
  double minimum_distance_change = 1.0e-4;
  /// \brief Increase factor for lambda if more damping needed.
  double lambda_increase = 5.0;
  /// \brief Decrease factor for lambda if new best solution is found.
  double lambda_decrease = 1.0 / 10.0;
  /// \brief Number of iterations to allow between calculating a new Jacobian.
  UnsignedIndex_t delay_jacobian_amount = 0;
  /// \brief Initial angle to use when first calculating Jacobian, equal to
  /// 5 degrees in radians.
  double initial_angle = 0.001 * 0.0174533;  // 1e-3 Deg in radians
  /// \brief Initial distance to use when first calculating Jacobian.
  double initial_distance = 0.001;
  /// \brief Angle change to use when calculating finite-difference Jacobian.
  double finite_difference_angle = 0.001 * 0.0174533;  // 1e-3 Deg in radians
};

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_OPTIMIZATION_BEHAVIOR_H_
