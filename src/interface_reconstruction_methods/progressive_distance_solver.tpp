// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_TPP_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_TPP_

#include <algorithm>

#include "src/generic_cutting/half_edge_cutting/half_edge_cutting.h"
#include "src/generic_cutting/half_edge_cutting/half_edge_cutting_initializer.tpp"

namespace IRL {

template <class CellType>
ProgressiveDistanceSolver<CellType>::ProgressiveDistanceSolver(
    const CellType& a_cell, const double a_volume_fraction,
    const double a_volume_fraction_tolerance,
    const PlanarSeparator& a_reconstruction)
    : target_volume_fraction_m(a_volume_fraction),
      volume_fraction_tolerance_m(a_volume_fraction_tolerance),
      reconstruction_m(a_reconstruction),
      flipped_solution_m(false) {
  assert(a_volume_fraction >= 0.0);
  assert(a_volume_fraction <= 1.0);
  assert(a_reconstruction.getNumberOfPlanes() == 1);
  this->solveForDistance(a_cell);
}

template <class CellType>
void ProgressiveDistanceSolver<CellType>::solve(
    const CellType& a_cell, const double a_volume_fraction,
    const double a_volume_fraction_tolerance,
    const PlanarSeparator& a_reconstruction) {
  assert(a_volume_fraction >= 0.0);
  assert(a_volume_fraction <= 1.0);
  assert(a_reconstruction.getNumberOfPlanes() == 1);
  target_volume_fraction_m = a_volume_fraction;
  volume_fraction_tolerance_m = a_volume_fraction_tolerance;
  reconstruction_m = a_reconstruction;
  flipped_solution_m = false;
  this->solveForDistance(a_cell);
}

template <class CellType>
SmallVector<double, 1> ProgressiveDistanceSolver<CellType>::getDistances(void) {
  return SmallVector<double, 1>({distances_m});
}

template <class CellType>
double ProgressiveDistanceSolver<CellType>::getDistances(
    const UnsignedIndex_t a_p) {
  assert(a_p == 0);
  return distances_m;
}

template <class CellType>
void ProgressiveDistanceSolver<CellType>::setup(const CellType& a_cell) {
  if (target_volume_fraction_m > 0.5) {
    flipped_solution_m = true;
    target_volume_fraction_m = 1.0 - target_volume_fraction_m;
    reconstruction_m[0] = reconstruction_m[0].generateFlippedPlane();
  }

  initial_cell_volume_m = a_cell.calculateVolume();

  const double max_dist =
      reconstruction_m[0].normal() * a_cell.calculateCentroid();
  UnsignedIndex_t accepted_dist = 0;
  sorted_distances_m.resize(a_cell.getNumberOfVertices());
  for (const auto& vertex : a_cell) {
    const double dist = reconstruction_m[0].normal() * vertex;
    if (dist <= max_dist) {
      sorted_distances_m[accepted_dist] = dist;
      ++accepted_dist;
    }
  }
  sorted_distances_m[accepted_dist] = max_dist;
  sorted_distances_m.resize(accepted_dist + 1);
  // Sort into ascending
  std::sort(sorted_distances_m.begin(), sorted_distances_m.end());

  assert(this->isBoundsTrueBounds(a_cell));
}

template <class CellType>
bool ProgressiveDistanceSolver<CellType>::isBoundsTrueBounds(
    const CellType& a_cell) {
  for (const auto& plane : reconstruction_m) {
    // Make sure plane has a valid normal
    assert(magnitude(plane.normal()) > 0.9);
  }
  reconstruction_m[0].distance() = sorted_distances_m.front();
  double result = std::fabs(
      getVolumeFraction<RecursiveSimplexCutting>(a_cell, reconstruction_m));
  if (result > target_volume_fraction_m) {
    // Volume fraction should be less than target_volume_fraction_m
    std::cout << "Lower bound result failed: " << result << std::endl;
    return false;
  }
  reconstruction_m[0].distance() = sorted_distances_m.back();
  result = target_volume_fraction_m -
           getVolumeFraction<RecursiveSimplexCutting>(a_cell, reconstruction_m);
  if (result > 0.0) {
    // Volume fraction should be > target_volume_fraction_m.
    std::cout << "Upper bound result failed: " << result << std::endl;
    return false;
  }
  return true;
}

template <class CellType>
void ProgressiveDistanceSolver<CellType>::solveForDistance(
    const CellType& a_cell) {
  // Setup the system
  this->setup(a_cell);

  // Generate half-edge structure
  // Reuse memory used in generic cutting.
  auto& complete_polytope = setHalfEdgeStructure(a_cell);
  auto under_polytope = generateSegmentedVersion<CellType>(&complete_polytope);
  decltype(under_polytope) over_polytope;

  // Keep bisectioning between nodes until we find between where the solution
  // lays. Will then just have a prismatoid that needs to be optimized over
  std::array<UnsignedIndex_t, 3> bounding_indices{
      {0, 0, static_cast<UnsignedIndex_t>(sorted_distances_m.size()) - 1}};
  std::array<double, 3> bounding_values{{0.0, 0.0, 0.5}};
  double found_vof_amount = 0.0;
  while (bounding_indices[2] - bounding_indices[1] > 1) {
    bounding_indices[1] = (bounding_indices[0] + bounding_indices[2]) / 2;
    reconstruction_m[0].distance() = sorted_distances_m[bounding_indices[1]];
    splitHalfEdgePolytope(&under_polytope, &over_polytope, &complete_polytope,
                          reconstruction_m[0]);
    const double volume_fraction_cut =
        under_polytope.calculateVolume() / initial_cell_volume_m;
    const double volume_fraction_guess = volume_fraction_cut + found_vof_amount;
    if (volume_fraction_guess >
        target_volume_fraction_m + volume_fraction_tolerance_m) {
      bounding_indices[2] = bounding_indices[1];
      bounding_values[2] = volume_fraction_guess;
    } else if (volume_fraction_guess <
               target_volume_fraction_m - volume_fraction_tolerance_m) {
      bounding_indices[0] = bounding_indices[1];
      bounding_values[0] = volume_fraction_guess;
      under_polytope = over_polytope;
      found_vof_amount += volume_fraction_cut;
    } else {
      distances_m = reconstruction_m[0].distance();
      distances_m *= flipped_solution_m ? -1.0 : 1.0;
      return;
    }
  }

  // Perform Brent's method  on remaining prismatoid
  // double a_delta = 1.0e-8;
  // double a_bracket_0 = sorted_distances_m[bounding_indices[0]];
  // double a_bracket_1 = sorted_distances_m[bounding_indices[2]];
  // double error_bracket_0 = bounding_values[0] - target_volume_fraction_m;
  // double error_bracket_1 = bounding_values[2] - target_volume_fraction_m;
  // if (std::fabs(error_bracket_0) < std::fabs(error_bracket_1)) {
  //   std::swap(a_bracket_0, a_bracket_1);
  //   std::swap(error_bracket_0, error_bracket_1);
  // }
  // double c = a_bracket_0;
  // double error_c = error_bracket_0;
  // double d = 0.0;
  // double s = 0.0;
  // double error_s = error_bracket_1;
  // bool mflag = true;
  // double error = std::fabs(error_s);
  // while (error > volume_fraction_tolerance_m) {
  //   if (error_bracket_0 != error_c && error_bracket_1 != error_c) {
  //     s = a_bracket_0 * error_bracket_1 * error_c /
  //             ((error_bracket_0 - error_bracket_1) *
  //              (error_bracket_0 - error_c)) +
  //         a_bracket_1 * error_bracket_0 * error_c /
  //             ((error_bracket_1 - error_bracket_0) *
  //              (error_bracket_1 - error_c)) +
  //         c * error_bracket_0 * error_bracket_1 /
  //             ((error_c - error_bracket_0) * (error_c - error_bracket_1));
  //   } else {
  //     s = a_bracket_1 - error_bracket_1 * (a_bracket_1 - a_bracket_0) /
  //                           (error_bracket_1 - error_bracket_0);
  //   }
  //   if (s < 0.25 * (3.0 * a_bracket_0 + a_bracket_1) || s > a_bracket_1 ||
  //       (mflag &&
  //        std::fabs(s - a_bracket_1) >= 0.5 * std::fabs(a_bracket_1 - c)) ||
  //       (!mflag && std::fabs(s - a_bracket_1) >= 0.5 * std::fabs(c - d)) ||
  //       (mflag && std::fabs(a_bracket_1 - c) < a_delta) ||
  //       (!mflag && std::fabs(c - d) < a_delta)) {
  //     s = 0.5 * (a_bracket_0 + a_bracket_1);
  //     mflag = true;
  //   } else {
  //     mflag = false;
  //   }
  //   reconstruction_m[0].distance() = s;
  //   splitHalfEdgePolytope(&under_polytope, &over_polytope,
  //   &complete_polytope,
  //                         reconstruction_m[0]);
  //   const double volume_fraction_cut =
  //       under_polytope.calculateVolume() / initial_cell_volume_m;
  //   const double volume_fraction_guess = volume_fraction_cut +
  //   found_vof_amount; error_s = volume_fraction_guess -
  //   target_volume_fraction_m; if (error_s < 0.0) {
  //     under_polytope = over_polytope;
  //     found_vof_amount += volume_fraction_cut;
  //   }
  //   d = c;
  //   c = a_bracket_1;
  //   error_c = error_bracket_1;
  //   if (error_bracket_0 * error_s < 0.0) {
  //     a_bracket_1 = s;
  //     error_bracket_1 = error_s;
  //   } else {
  //     a_bracket_0 = s;
  //     error_bracket_0 = error_s;
  //   }
  //   if (std::fabs(error_bracket_0) < std::fabs(error_bracket_1)) {
  //     std::swap(a_bracket_0, a_bracket_1);
  //     std::swap(error_bracket_0, error_bracket_1);
  //   }
  //   error = {std::min(std::fabs(error_bracket_1), std::fabs(error_s))};
  // }
  // distances_m =
  //     std::fabs(error_s) < std::fabs(error_bracket_1) ? s : a_bracket_1;
  // distances_m *= flipped_solution_m ? -1.0 : 1.0;
  // return;

  // // Now just have a prismatoid left. Perform secant, fall back to
  // Bisection if it fails.
  reconstruction_m[0].distance() = sorted_distances_m[bounding_indices[0]];
  double delta = sorted_distances_m[bounding_indices[0]] -
                 sorted_distances_m[bounding_indices[2]];
  double old_error = bounding_values[2] - target_volume_fraction_m;
  double error = bounding_values[0] - target_volume_fraction_m;

  for (UnsignedIndex_t iter = 0; iter < max_iter_m; ++iter) {
    delta *= -error / safelyEpsilon(error - old_error);
    old_error = error;
    reconstruction_m[0].distance() += delta;
    splitHalfEdgePolytope(&under_polytope, &over_polytope, &complete_polytope,
                          reconstruction_m[0]);
    const double volume_fraction_cut =
        under_polytope.calculateVolume() / initial_cell_volume_m;
    const double volume_fraction_guess = volume_fraction_cut + found_vof_amount;
    error = volume_fraction_guess - target_volume_fraction_m;

    if (volume_fraction_guess >
        target_volume_fraction_m + volume_fraction_tolerance_m) {
      if (reconstruction_m[0].distance() < bounding_values[2]) {
        bounding_values[2] = reconstruction_m[0].distance();
      }
    } else if (volume_fraction_guess <
               target_volume_fraction_m - volume_fraction_tolerance_m) {
      if (reconstruction_m[0].distance() > bounding_values[0]) {
        bounding_values[0] = reconstruction_m[0].distance();
      }
      under_polytope = over_polytope;
      found_vof_amount += volume_fraction_cut;
    } else {
      distances_m = reconstruction_m[0].distance();
      distances_m *= flipped_solution_m ? -1.0 : 1.0;
      return;
    }
  }

  // Perform bisection since secant failed to find answer within tolerance.
  for (UnsignedIndex_t iter = 0; iter < max_bisection_iter; ++iter) {
    bounding_values[1] = 0.5 * (bounding_values[0] + bounding_values[2]);
    reconstruction_m[0].distance() = bounding_values[1];
    splitHalfEdgePolytope(&under_polytope, &over_polytope, &complete_polytope,
                          reconstruction_m[0]);
    const double volume_fraction_cut =
        under_polytope.calculateVolume() / initial_cell_volume_m;
    const double volume_fraction_guess = volume_fraction_cut + found_vof_amount;
    if (volume_fraction_guess >
        target_volume_fraction_m + volume_fraction_tolerance_m) {
      bounding_values[2] = bounding_values[1];
    } else if (volume_fraction_guess <
               target_volume_fraction_m - volume_fraction_tolerance_m) {
      bounding_values[0] = bounding_values[1];
      under_polytope = over_polytope;
      found_vof_amount += volume_fraction_cut;
    } else {
      distances_m = reconstruction_m[0].distance();
      distances_m *= flipped_solution_m ? -1.0 : 1.0;
      return;
    }
  }
}  // namespace IRL

}  // namespace IRL

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_TPP_
