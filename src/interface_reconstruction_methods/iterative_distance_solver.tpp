// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_ITERATIVE_DISTANCE_SOLVER_TPP_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_ITERATIVE_DISTANCE_SOLVER_TPP_

#include "src/optimization/brents_method.h"

namespace IRL {

template <class CellType, UnsignedIndex_t kMaxPlanes>
IterativeSolverForDistance<CellType, kMaxPlanes>::IterativeSolverForDistance(
    const CellType& a_cell, const double a_volume_fraction,
    const double a_volume_fraction_tolerance,
    const PlanarSeparator& a_reconstruction)
    : cell_m(&a_cell),
      target_volume_fraction_m(a_volume_fraction),
      volume_fraction_tolerance_m(a_volume_fraction_tolerance),
      reconstruction_m(a_reconstruction),
      flipped_solution_m(false) {
  this->solveForDistance();
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CellType, kMaxPlanes>::solve(
    const CellType& a_cell, const double a_volume_fraction,
    const double a_volume_fraction_tolerance,
    const PlanarSeparator& a_reconstruction) {
  cell_m = &a_cell;
  target_volume_fraction_m = a_volume_fraction;
  volume_fraction_tolerance_m = a_volume_fraction_tolerance;
  reconstruction_m = a_reconstruction;
  flipped_solution_m = false;
  this->solveForDistance();
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CellType, kMaxPlanes>::updateGuess(
    double* a_delta) {
  current_guess_m += *a_delta;
  for (UnsignedIndex_t n = 0; n < distances_m.size(); ++n) {
    distances_m[n] = initial_distances_m[n] + current_guess_m;
  }
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CellType, kMaxPlanes>::setGuess(
    double* a_guess) {
  current_guess_m = *a_guess;
  for (UnsignedIndex_t n = 0; n < distances_m.size(); ++n) {
    distances_m[n] = initial_distances_m[n] + current_guess_m;
  }
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
double
IterativeSolverForDistance<CellType, kMaxPlanes>::calculateSignedScalarError(
    void) {
  reconstruction_m.setDistances(distances_m);
  return target_volume_fraction_m -
         getVolumeFraction<RecursiveSimplexCutting>(*cell_m, reconstruction_m);
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
bool IterativeSolverForDistance<CellType, kMaxPlanes>::errorTooHigh(
    const double a_error) {
  return std::fabs(a_error) > volume_fraction_tolerance_m;
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
bool IterativeSolverForDistance<CellType, kMaxPlanes>::iterationTooHigh(
    const UnsignedIndex_t a_iteration) {
  return a_iteration > max_iter_m;
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
const SmallVector<double, kMaxPlanes>&
IterativeSolverForDistance<CellType, kMaxPlanes>::getDistances(void) {
  return distances_m;
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
double IterativeSolverForDistance<CellType, kMaxPlanes>::getDistances(
    const UnsignedIndex_t a_p) {
  assert(a_p < distances_m.size());
  return distances_m[a_p];
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CellType, kMaxPlanes>::setup(void) {
  if (target_volume_fraction_m > 0.5) {
    flipped_solution_m = true;
    target_volume_fraction_m = 1.0 - target_volume_fraction_m;
    for (auto& plane : reconstruction_m) {
      plane = plane.generateFlippedPlane();
    }
    if (reconstruction_m.getNumberOfPlanes() > 1) {
      double new_flip = reconstruction_m.flip() == 1.0 ? -1.0 : 1.0;
      reconstruction_m.setFlip(new_flip);
    }
  }

  initial_distances_m.resize(reconstruction_m.getNumberOfPlanes());
  distances_m.resize(reconstruction_m.getNumberOfPlanes());
  this->checkIfStaticAllocationExceeded();
  for (UnsignedIndex_t n = 0; n < reconstruction_m.getNumberOfPlanes(); ++n) {
    initial_distances_m[n] = reconstruction_m[n].distance();
  }
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
bool IterativeSolverForDistance<CellType, kMaxPlanes>::isBoundsTrueBounds(
    void) {
  for (const auto& plane : reconstruction_m) {
    // Make sure plane has a valid normal
    assert(magnitude(plane.normal()) > 0.9);
  }
  double guess[1];
  guess[0] = bound_value_m[0];
  this->setGuess(guess);
  reconstruction_m.setDistances(distances_m);
  double result = std::fabs(
      getVolumeFraction<RecursiveSimplexCutting>(*cell_m, reconstruction_m));
  if (result > 1.0e-14) {
    // Volume fraction should be 0.0.
    std::cout << "Lower bound result failed: " << result << std::endl;
    return false;
  }
  guess[0] = bound_value_m[2];
  this->setGuess(guess);
  reconstruction_m.setDistances(distances_m);
  result = std::fabs(1.0 - getVolumeFraction<RecursiveSimplexCutting>(
                               *cell_m, reconstruction_m));
  if (result > 1.0e-14) {
    // Volume fraction should be 1.0.
    std::cout << "Upper bound result failed: " << result << std::endl;
    return false;
  }
  guess[0] = 0.0;
  this->setGuess(guess);
  return true;
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CellType, kMaxPlanes>::calculateBounds(void) {
  double lower_bound = DBL_MAX;
  double upper_bound = -DBL_MAX;
  for (const auto& plane : reconstruction_m) {
    for (const auto& vertex : (*cell_m)) {
      lower_bound = std::min(lower_bound, vertex * plane.normal());
      upper_bound = std::max(upper_bound, vertex * plane.normal());
    }
  }
  lower_bound = lower_bound - 0.1 * std::fabs(lower_bound);
  upper_bound = upper_bound + 0.1 * std::fabs(upper_bound);

  bound_value_m[0] = DBL_MAX;
  bound_value_m[2] = -DBL_MAX;
  this->checkIfStaticAllocationExceeded();
  for (UnsignedIndex_t n = 0; n < reconstruction_m.getNumberOfPlanes(); ++n) {
    bound_value_m[0] =
        std::min(bound_value_m[0], lower_bound - initial_distances_m[n]);
    bound_value_m[2] =
        std::max(bound_value_m[2], upper_bound - initial_distances_m[n]);
  }
  bound_value_m[1] = 0.0;

  bound_error_m[0] = target_volume_fraction_m;
  bound_error_m[1] = 0.0;
  bound_error_m[2] = target_volume_fraction_m - 1.0;

  assert(isBoundsTrueBounds());
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CellType, kMaxPlanes>::solveForDistance(void) {
  // Setup the system
  this->setup();

  // Try to solve using secant method
  Secant<IterativeSolverForDistance> secant_solver;
  secant_solver.solve(this,
                      0.01 * std::pow(cell_m->calculateVolume(), 1.0 / 3.0));
  if (secant_solver.getReason() >= 0) {
    if (flipped_solution_m) {
      for (auto& member : distances_m) {
        member = -member;
      }
    }
    return;
  }

  //  If made it this far, Secant method wasn't good enough
  //  Now lets default to bisection to finish this out
  this->calculateBounds();
  Bisection<IterativeSolverForDistance> bisection_solver;
  bisection_solver.solve(this, bound_value_m[0], bound_value_m[2],
                         max_bisection_iter);
  if (flipped_solution_m) {
    for (auto& member : distances_m) {
      member = -member;
    }
  }
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CellType, kMaxPlanes>::updateBoundaries(
    const double a_error) {
  // Check for better negative bracket end
  if (a_error > 0.0 && a_error < bound_error_m[0] &&
      current_guess_m > bound_value_m[0]) {
    bound_value_m[0] = current_guess_m;
    bound_error_m[0] = a_error;
  }
  // Check for better positive bracket end
  if (a_error < 0.0 && a_error > bound_error_m[2] &&
      current_guess_m < bound_value_m[2]) {
    bound_value_m[2] = current_guess_m;
    bound_error_m[2] = a_error;
  }
}

template <class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<
    CellType, kMaxPlanes>::checkIfStaticAllocationExceeded(void) const {
#ifndef NDEBUG_PERF
  if ((initial_distances_m.capacity() > kMaxPlanes) ||
      (distances_m.capacity() > kMaxPlanes)) {
    std::cout << "Static allocation size for SmallVector exceeded in "
                 "IterativeSolverForDistance. Expect performance "
                 "penalty if this happens frequently."
              << std::endl;
  }
#endif
}

}  // namespace IRL

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_ITERATIVE_DISTANCE_SOLVER_TPP_
