// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_PARABOLOID_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_PARABOLOID_H_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/generic_cutting/paraboloid_intersection/moment_contributions.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/helpers/helper.h"
#include "irl/helpers/mymath.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/parameters/constants.h"

namespace IRL {

/// \brief Volume conserving distance-finding routine for
/// single paraboloid reconstructions.
///
/// \param[in] volume_fraction_m Volume fraction to recreate
/// \param[in] volume_fraction_tolerance_m Tolerance to recreate
/// `a_volume_fraction` within
/// \param[in] reconstruction_m A copy of the reconstruction
/// to find distances for.
/// \param[out] distances_m Correct distance to plane is stored
/// in distances_m after construction of the object.
template <class CellType>
class ProgressiveDistanceSolverParaboloid {
  /// \brief Max number of iterations for the secant Solver
  static constexpr UnsignedIndex_t max_iter_m = {80};
  static constexpr UnsignedIndex_t max_bisection_iter = {80};

 public:
  /// \brief Constructor that initializes the class for optimization
  /// and solves for the distance.

  /// \brief Default constructor
  ProgressiveDistanceSolverParaboloid(void) = default;

  /// \brief Constructor that loads in the neccessary information to
  /// the class and then solves for volume-conserving distance to two planes.
  ///
  /// \param[in] a_volume_fraction Volume fraction to match.
  /// \param[in] a_volume_fraction_tolerance Tolerance to recreate volume
  /// fraction within.
  /// \param[in] a_reconstruction Reconstruction to find distances for.
  ProgressiveDistanceSolverParaboloid(const CellType& a_cell,
                                      const double a_volume_fraction,
                                      const double a_volume_fraction_tolerance,
                                      const Paraboloid& a_reconstruction);

  /// \brief Reinitialize solver and solve for distance
  /// to each plane.
  ///
  /// \param[in] a_volume_fraction Volume fraction to match.
  /// \param[in] a_volume_fraction_tolerance Tolerance to recreate volume
  /// fraction within.
  /// \param[in] a_reconstruction Reconstruction to find distances for.
  void solve(const CellType& a_cell, const double a_volume_fraction,
             const double a_volume_fraction_tolerance,
             const Paraboloid& a_reconstruction);

  /// \brief Return pointer to distances_m to be used for changing a
  /// reconstruction.
  double getDistance(void);

  /// \brief Make sure initial bounds are set up correctly.
  inline bool isBoundsTrueBounds(const CellType& a_cell);

  /// \brief Default destructor
  ~ProgressiveDistanceSolverParaboloid(void) = default;

 private:
  /// \brief Main driving function that solves for distances
  inline void solveForDistance(const CellType& a_cell);

  /// \brief Starting volume of the cell
  double initial_cell_volume_m;
  /// \brief Target volume fraction to match
  double target_volume_fraction_m;
  /// \brief Tolerance to match volume fraction within
  double volume_fraction_tolerance_m;
  // \brief Distance to each vertex in cell
  std::vector<double> sorted_distances_m;
  /// \brief The reconstruction the distance is being found for
  Paraboloid reconstruction_m;
  /// \brief Distances that satisfy tolerance.
  double distances_m;
};

template <class SegmentedHalfEdgePolytopeType>
void shiftSegmentedPolytopeAndUpdateFaces(
    SegmentedHalfEdgePolytopeType* a_seg_half_edge, const double a_shift);

}  // namespace IRL

#include "irl/interface_reconstruction_methods/progressive_distance_solver_paraboloid.tpp"

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_PARABOLOID_H_
