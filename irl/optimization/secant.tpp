// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_OPTIMIZATION_SECANT_TPP_
#define IRL_OPTIMIZATION_SECANT_TPP_

namespace IRL {

template <class OptimizingClass>
void Secant<OptimizingClass>::solve(OptimizingClass* a_setup_otype,
                                    const double a_initial_delta) {
  otype_m = a_setup_otype;
  this->solve(a_initial_delta);
}

template <class OptimizingClass>
void Secant<OptimizingClass>::solve(const double a_initial_delta) {
  assert(otype_m != nullptr);

  delta_m = 0.0;
  otype_m->updateGuess(&delta_m);
  double old_error = {otype_m->calculateSignedScalarError()};

  // Take initial step and get new error
  delta_m = a_initial_delta;
  otype_m->updateGuess(&delta_m);
  double error = otype_m->calculateSignedScalarError();

  iteration_m = 0;
  while (otype_m->errorTooHigh(error)) {
    ++iteration_m;
    if (otype_m->iterationTooHigh(iteration_m)) {
      // Exiting because exceeding max iterations
      reason_for_exit_m = -1;
      return;
    }

    // Calculate update step
    delta_m = error * delta_m / safelyEpsilon(old_error - error);
    old_error = error;
    otype_m->updateGuess(&delta_m);
    error = otype_m->calculateSignedScalarError();
  }

  // Error driven below accetable level, set reason as iteration number
  reason_for_exit_m = static_cast<int>(iteration_m);
}

template <class OptimizingClass>
void Secant<OptimizingClass>::solveAndTrackBounds(
    OptimizingClass* a_setup_otype, const double a_initial_delta) {
  otype_m = a_setup_otype;
  this->solveAndTrackBounds(a_initial_delta);
}

template <class OptimizingClass>
void Secant<OptimizingClass>::solveAndTrackBounds(
    const double a_initial_delta) {
  assert(otype_m != nullptr);

  delta_m = 0.0;
  otype_m->updateGuess(&delta_m);
  double old_error = {otype_m->calculateSignedScalarError()};
  reason_for_exit_m = 0;

  // Take initial step and get new error
  delta_m = a_initial_delta;
  otype_m->updateGuess(&delta_m);
  double error = otype_m->calculateSignedScalarError();

  // Update boundaries that are being tracked
  otype_m->updateBoundaries(error);

  iteration_m = 0;
  while (otype_m->errorTooHigh(error)) {
    iteration_m++;
    if (otype_m->iterationTooHigh(iteration_m)) {
      // Exiting because exceeding max iterations
      reason_for_exit_m = -1;
      return;
    }

    // Calculate update step
    delta_m = error * delta_m / safelyEpsilon(old_error - error);
    old_error = error;
    otype_m->updateGuess(&delta_m);
    error = otype_m->calculateSignedScalarError();
    otype_m->updateBoundaries(error);
  }

  // Error driven below accetable level, set reason as iteration number
  reason_for_exit_m = static_cast<int>(iteration_m);
}

template <class OptimizingClass>
int Secant<OptimizingClass>::getReason(void) const {
  return reason_for_exit_m;
}

}  // namespace IRL

#endif // IRL_OPTIMIZATION_SECANT_TPP_
