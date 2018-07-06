// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_ITERATIVE_DISTANCE_SOLVER_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_ITERATIVE_DISTANCE_SOLVER_H_

#include <algorithm>
#include <cmath>
#include <iostream>

#include "src/data_structures/stack_vector.h"
#include "src/generic_cutting/generic_cutting.h"
#include "src/geometry/general/normal.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/helpers/helper.h"
#include "src/optimization/bisection.h"
#include "src/optimization/secant.h"
#include "src/parameters/constants.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

/// \brief Volume conserving distance-finding routine
/// for two-plane reconstructions wrapped in a class.
///
/// This class finds the distance to the two planes
/// provided in a_reconstruction to recreate `a_volume_fraction`
/// within the tolerance `a_find_volume_tolerance`. This is
/// first attempted by a Newton-Raphson optimization. If
/// the optimization does not reach within the tolerance
/// in a given number of iterations, it then defaults to
/// using bisection with the boundaries being the
/// current best solutions on both sides of the
/// goal volume_fraction.
///
/// \param[in] volume_fraction_m Volume fraction to recreate
/// \param[in] volume_fraction_tolerance_m Tolerance to recreate
/// `a_volume_fraction` within
/// \param[in] reconstruction_m A copy of the reconstruction
/// to find distances for.
/// param[out] distances_m Correct distance to plane is stored
/// in distances_m after construction of the object.
template <class CellType, UnsignedIndex_t kMaxPlanes =
                              global_constants::MAX_PLANAR_SEPARATOR_PLANES>
class IterativeSolverForDistance {
  /// \brief Max number of iterations for the Newton-Raphson Solver
  static constexpr UnsignedIndex_t max_iter_m = {15};
  static constexpr UnsignedIndex_t max_bisection_iter = {40};

 public:
  /// \brief Constructor that initializes the class for optimization
  /// and solves for the distance.

  /// \brief Default constructor
  IterativeSolverForDistance(void) = default;

  /// \brief Constructor that loads in the neccessary information to
  /// the class and then solves for volume-conserving distance to two planes.
  ///
  /// \param[in] a_volume_fraction Volume fraction to match.
  /// \param[in] a_volume_fraction_tolerance Tolerance to recreate volume
  /// fraction within.
  /// \param[in] a_reconstruction Reconstruction to find distances for.
  IterativeSolverForDistance(const CellType& a_cell,
                             const double a_volume_fraction,
                             const double a_volume_fraction_tolerance,
                             const PlanarSeparator& a_reconstruction);
  /// \brief Reinitialize solver and solve for distance
  /// to each plane.
  ///
  /// \param[in] a_volume_fraction Volume fraction to match.
  /// \param[in] a_volume_fraction_tolerance Tolerance to recreate volume
  /// fraction within.
  /// \param[in] a_reconstruction Reconstruction to find distances for.
  void solve(const CellType& a_cell, const double a_volume_fraction,
             const double a_volume_fraction_tolerance,
             const PlanarSeparator& a_reconstruction);

  /// \brief Update the current guess and solution
  /// given the delta change.
  void updateGuess(double* a_delta);

  /// \brief  Set new guess distances to each plane.
  void setGuess(double* a_guess);

  /// \brief Calculate current volume fraction for reconstruction_m
  /// with the plane distances_m, return error
  double calculateSignedScalarError(void);

  /// \brief Return bool for whether the error is still too high.
  bool errorTooHigh(const double a_error);

  /// \brief Return bool for whether maximum iterations has been exceeded
  bool iterationTooHigh(const UnsignedIndex_t a_iteration);

  /// \brief Update known boundaries for solution.
  inline void updateBoundaries(const double a_error);

  /// \brief Return pointer to distances_m to be used for changing a
  /// reconstruction.
  const SmallVector<double, kMaxPlanes>& getDistances(void);

  /// \brief Return the distance for plane `a_p`.
  double getDistances(const UnsignedIndex_t a_p);

  /// \brief Default destructor
  ~IterativeSolverForDistance(void) = default;

 private:
  /// \brief Prepare the object to be passed to generalized solvers
  inline void setup(void);

  /// \brief Calculates max and min bounds for entire cell.
  inline void calculateBounds(void);

  /// \brief Main driving function that solves for distances
  inline void solveForDistance(void);

  /// \brief Make sure initial bounds are set up correctly.
  inline bool isBoundsTrueBounds(void);

  void checkIfStaticAllocationExceeded(void) const;

  /// \brief Cell distance is being calculated for.
  const CellType* cell_m;
  /// \brief Target volume fraction to match
  double target_volume_fraction_m;
  /// \brief Tolerance to match volume fraction within
  double volume_fraction_tolerance_m;
  /// \brief The reconstruction the distance is being found for
  PlanarSeparator reconstruction_m;
  /// \brief Clipped initial distance that was in reconstruction_m that the
  /// optimization will start from.
  SmallVector<double, kMaxPlanes> initial_distances_m;
  /// \brief Distances that the planes in reconstruction should be set to
  SmallVector<double, kMaxPlanes> distances_m;
  /// \brief The current solution
  double current_guess_m;
  bool flipped_solution_m;

  /// \brief Values bracketing solution
  std::array<double, 3> bound_value_m;
  /// \brief Errors for values bracketing solution
  std::array<double, 3> bound_error_m;
};

}  // namespace IRL

#include "src/interface_reconstruction_methods/iterative_distance_solver.tpp"

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_ITERATIVE_DISTANCE_SOLVER_H_
