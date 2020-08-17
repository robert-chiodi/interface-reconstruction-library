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
#include "src/optimization/illinois.h"

namespace IRL {

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
IterativeSolverForDistance<CuttingMethod, CellType, kMaxPlanes>::
    IterativeSolverForDistance(const CellType &a_cell,
                               const double a_volume_fraction,
                               const double a_volume_fraction_tolerance,
                               const PlanarSeparator &a_reconstruction)
    : cell_m(&a_cell), target_volume_fraction_m(a_volume_fraction),
      volume_fraction_tolerance_m(a_volume_fraction_tolerance),
      reconstruction_m(a_reconstruction), flipped_solution_m(false) {
  this->solveForDistance();
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CuttingMethod, CellType, kMaxPlanes>::solve(
    const CellType &a_cell, const double a_volume_fraction,
    const double a_volume_fraction_tolerance,
    const PlanarSeparator &a_reconstruction) {
  cell_m = &a_cell;
  target_volume_fraction_m = a_volume_fraction;
  volume_fraction_tolerance_m = a_volume_fraction_tolerance;
  reconstruction_m = a_reconstruction;
  flipped_solution_m = false;
  this->solveForDistance();
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CuttingMethod, CellType,
                                kMaxPlanes>::updateGuess(double *a_delta) {
  current_guess_m += *a_delta * characteristic_length_m;
  for (UnsignedIndex_t n = 0; n < distances_m.size(); ++n) {
    distances_m[n] = initial_distances_m[n] + current_guess_m;
  }
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CuttingMethod, CellType, kMaxPlanes>::setGuess(
    double *a_guess) {
  const double adjustment = *a_guess * characteristic_length_m;
  for (UnsignedIndex_t n = 0; n < distances_m.size(); ++n) {
    distances_m[n] = initial_distances_m[n] + adjustment;
  }
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
double
IterativeSolverForDistance<CuttingMethod, CellType,
                           kMaxPlanes>::calculateSignedScalarError(void) {
  reconstruction_m.setDistances(distances_m);
  return target_volume_fraction_m -
         getVolumeFraction<CuttingMethod>(*cell_m, reconstruction_m);
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
bool IterativeSolverForDistance<
    CuttingMethod, CellType, kMaxPlanes>::errorTooHigh(const double a_error) {
  return std::fabs(a_error) > volume_fraction_tolerance_m;
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
bool IterativeSolverForDistance<CuttingMethod, CellType, kMaxPlanes>::
    iterationTooHigh(const UnsignedIndex_t a_iteration) {
  return a_iteration > max_iter_m;
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
const SmallVector<double, kMaxPlanes> &
IterativeSolverForDistance<CuttingMethod, CellType, kMaxPlanes>::getDistances(
    void) const {
  return distances_m;
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
double
IterativeSolverForDistance<CuttingMethod, CellType, kMaxPlanes>::getDistances(
    const UnsignedIndex_t a_p) const {
  assert(a_p < distances_m.size());
  return distances_m[a_p];
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CuttingMethod, CellType, kMaxPlanes>::setup(
    void) {
  characteristic_length_m = std::cbrt(cell_m->calculateVolume());
  current_guess_m = 0.0;
  if (target_volume_fraction_m > 0.5) {
    flipped_solution_m = true;
    target_volume_fraction_m = 1.0 - target_volume_fraction_m;
    for (auto &plane : reconstruction_m) {
      plane = plane.generateFlippedPlane();
    }
    if(reconstruction_m.getNumberOfPlanes() > 1){
      reconstruction_m.setFlip(reconstruction_m.flip() == 1.0 ? -1.0 : 1.0);
    }
  }

  initial_distances_m.resize(reconstruction_m.getNumberOfPlanes());
  distances_m.resize(reconstruction_m.getNumberOfPlanes());
  this->checkIfStaticAllocationExceeded();
  for (UnsignedIndex_t n = 0; n < reconstruction_m.getNumberOfPlanes(); ++n) {
    initial_distances_m[n] = reconstruction_m[n].distance();
  }
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
bool IterativeSolverForDistance<CuttingMethod, CellType,
                                kMaxPlanes>::isBoundsTrueBounds(void) {
  for (const auto &plane : reconstruction_m) {
    // Make sure plane has a valid normal
    assert(magnitude(plane.normal()) > 0.9);
  }
  double guess[1];
  guess[0] = bound_value_m[0];
  this->setGuess(guess);
  reconstruction_m.setDistances(distances_m);
  double result =
      std::fabs(getVolumeFraction<CuttingMethod>(*cell_m, reconstruction_m));
  if (result > 1.0e-14) {
    // Volume fraction should be 0.0.
    std::cout << "Lower bound result failed: " << result << std::endl;
    return false;
  }
  guess[0] = bound_value_m[1];
  this->setGuess(guess);
  reconstruction_m.setDistances(distances_m);
  result = std::fabs(
      1.0 - getVolumeFraction<CuttingMethod>(*cell_m, reconstruction_m));
  if (result > 1.0e-14) {
    // Volume fraction should be 1.0.
    std::cout << "Upper bound result failed: " << result << std::endl;
    return false;
  }
  guess[0] = 0.0;
  this->setGuess(guess);
  return true;
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CuttingMethod, CellType,
                                kMaxPlanes>::calculateBounds(void) {
  double lower_bound = DBL_MAX;
  double upper_bound = -DBL_MAX;
  for (const auto &plane : reconstruction_m) {
    for (const auto &vertex : (*cell_m)) {
      lower_bound = std::min(lower_bound, vertex * plane.normal());
      upper_bound = std::max(upper_bound, vertex * plane.normal());
    }
  }
  lower_bound -= 1.0e-8 * std::fabs(lower_bound);
  upper_bound += 1.0e-8 * std::fabs(upper_bound);

  for (auto &distance : initial_distances_m) {
    distance = clipBetween(lower_bound, distance, upper_bound);
  }

  bound_value_m[0] = DBL_MAX;
  bound_value_m[1] = -DBL_MAX;
  this->checkIfStaticAllocationExceeded();
  for (UnsignedIndex_t n = 0; n < reconstruction_m.getNumberOfPlanes(); ++n) {
    bound_value_m[0] =
        std::min(bound_value_m[0], lower_bound - initial_distances_m[n]);
    bound_value_m[1] =
        std::max(bound_value_m[1], upper_bound - initial_distances_m[n]);
  }
  bound_value_m[0] /= characteristic_length_m;
  bound_value_m[1] /= characteristic_length_m;
  assert(isBoundsTrueBounds());
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CuttingMethod, CellType,
                                kMaxPlanes>::solveForDistance(void) {
  // Setup the system
  this->setup();

  this->calculateBounds();
  Illinois<IterativeSolverForDistance> illinois_solver;
  illinois_solver.solve(this, bound_value_m[0], bound_value_m[1]);
  if (flipped_solution_m) {
    for (auto &member : distances_m) {
      member = -member;
    }
  }
}

template <class CuttingMethod, class CellType, UnsignedIndex_t kMaxPlanes>
void IterativeSolverForDistance<CuttingMethod, CellType, kMaxPlanes>::
    checkIfStaticAllocationExceeded(void) const {
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

} // namespace IRL

#endif // SRC_INTERFACE_RECONSTRUCTION_METHODS_ITERATIVE_DISTANCE_SOLVER_TPP_
