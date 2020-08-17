// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_H_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "src/data_structures/small_vector.h"
#include "src/generic_cutting/generic_cutting.h"
#include "src/geometry/general/normal.h"
#include "src/geometry/general/plane.h"
#include "src/helpers/helper.h"
#include "src/parameters/constants.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

/// \brief Volume conserving distance-finding routine for
/// single plane reconstructions.
///
/// \param[in] volume_fraction_m Volume fraction to recreate
/// \param[in] volume_fraction_tolerance_m Tolerance to recreate
/// `a_volume_fraction` within
/// \param[in] reconstruction_m A copy of the reconstruction
/// to find distances for.
/// param[out] distances_m Correct distance to plane is stored
/// in distances_m after construction of the object.
template <class CellType>
class ProgressiveDistanceSolver {
  /// \brief Max number of iterations for the Newton-Raphson Solver
  static constexpr UnsignedIndex_t max_iter_m = {15};
  static constexpr UnsignedIndex_t max_bisection_iter = {80};

 public:
  /// \brief Constructor that initializes the class for optimization
  /// and solves for the distance.

  /// \brief Default constructor
  ProgressiveDistanceSolver(void) = default;

  /// \brief Constructor that loads in the neccessary information to
  /// the class and then solves for volume-conserving distance to two planes.
  ///
  /// \param[in] a_volume_fraction Volume fraction to match.
  /// \param[in] a_volume_fraction_tolerance Tolerance to recreate volume
  /// fraction within.
  /// \param[in] a_reconstruction Reconstruction to find distances for.
  ProgressiveDistanceSolver(const CellType& a_cell,
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

  /// \brief Return pointer to distances_m to be used for changing a
  /// reconstruction.
  SmallVector<double, 1> getDistances(void);

  /// \brief Return the distance for plane `a_p`.
  double getDistances(const UnsignedIndex_t a_p);

  /// \brief Default destructor
  ~ProgressiveDistanceSolver(void) = default;

 private:
  /// \brief Prepare the object to be passed to generalized solvers
  inline void setup(const CellType& a_cell);

  /// \brief Main driving function that solves for distances
  inline void solveForDistance(const CellType& a_cell);

  /// \brief Make sure initial bounds are set up correctly.
  inline bool isBoundsTrueBounds(const CellType& a_cell);

  /// \brief Starting volume of the cell
  double initial_cell_volume_m;
  /// \brief Target volume fraction to match
  double target_volume_fraction_m;
  /// \brief Tolerance to match volume fraction within
  double volume_fraction_tolerance_m;
  /// \brief The reconstruction the distance is being found for
  PlanarSeparator reconstruction_m;
  /// \brief Distances that satisfy tolerance.
  double distances_m;
  // \brief Distance to each vertex in cell
  std::vector<double> sorted_distances_m;
  /// \brief Bool indicating whether solution was found for flipped phase.
  bool flipped_solution_m;
};

}  // namespace IRL

#include "src/interface_reconstruction_methods/progressive_distance_solver.tpp"

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_H_
