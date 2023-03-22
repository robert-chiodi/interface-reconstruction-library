// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_CONSTRAINED_OPTIMIZATION_BEHAVIOR_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_CONSTRAINED_OPTIMIZATION_BEHAVIOR_H_

#include "irl/parameters/defined_types.h"

namespace IRL {

/// \brief Struct to contain PLVIRA/PMoF optimization parameters
struct ConstrainedOptimizationBehavior {
  /// \brief If `this->calculateScalarError()` is less than this, exit.
  double acceptable_squared_error = 1.0e-12 * 1.0e-12;
  /// \brief If `this->calculateScalarError()` is less than this, exit.
  double acceptable_max_error = 1.0e-8;
  /// \brief If `this->calculateScalarError()` is less than this, exit.
  double acceptable_squared_constraint_error = 1.0e-5 * 1.0e-5;
  /// \brief Maximum number of attempted iterations before exiting.
  UnsignedIndex_t maximum_iterations = 20;
  /// \brief Maximum number of attempted sub-iterations before exiting.
  UnsignedIndex_t maximum_sub_iterations = 25;
  /// \brief Factor for increasing the penalty parameter of the LM algo.
  double penalty_param_increase_factor = 1.8;
  /// \brief Factor for increasing the damping parameter of the LM algo.
  double damping_param_increase_factor = 5.0;
  /// \brief Initial factor for the damping parameter of the LM algo.
  double damping_param_initial_factor = 1.0;
};

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_CONSTRAINED_OPTIMIZATION_BEHAVIOR_H_
