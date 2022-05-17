// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_PARABOLOID_TPP_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_PARABOLOID_TPP_

#include <algorithm>

#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_initializer.tpp"

namespace IRL {

template <class CellType>
ProgressiveDistanceSolverParaboloid<CellType>::
    ProgressiveDistanceSolverParaboloid(
        const CellType& a_cell, const double a_volume_fraction,
        const double a_volume_fraction_tolerance,
        const Paraboloid& a_reconstruction)
    : target_volume_fraction_m(a_volume_fraction),
      volume_fraction_tolerance_m(a_volume_fraction_tolerance),
      reconstruction_m(a_reconstruction) {
  assert(a_volume_fraction >= 0.0);
  assert(a_volume_fraction <= 1.0);
  this->solveForDistance(a_cell);
}

template <class CellType>
void ProgressiveDistanceSolverParaboloid<CellType>::solve(
    const CellType& a_cell, const double a_volume_fraction,
    const double a_volume_fraction_tolerance,
    const Paraboloid& a_reconstruction) {
  assert(a_volume_fraction >= 0.0);
  assert(a_volume_fraction <= 1.0);
  target_volume_fraction_m = a_volume_fraction;
  volume_fraction_tolerance_m = a_volume_fraction_tolerance;
  reconstruction_m = a_reconstruction;
  this->solveForDistance(a_cell);
}

template <class CellType>
double ProgressiveDistanceSolverParaboloid<CellType>::getDistance(void) {
  return distances_m;
}

template <class CellType>
bool ProgressiveDistanceSolverParaboloid<CellType>::isBoundsTrueBounds(
    const CellType& a_cell) {
  const auto datum = reconstruction_m.getDatum();
  const auto& ref_frame = reconstruction_m.getReferenceFrame();
  const auto& aligned_paraboloid = reconstruction_m.getAlignedParaboloid();

  const Pt low_datum = datum + sorted_distances_m.front() * ref_frame[2];
  std::cout << "Low datum = " << low_datum << std::endl;
  const Pt high_datum = datum + sorted_distances_m.back() * ref_frame[2];
  std::cout << "High datum = " << high_datum << std::endl;
  const auto low_paraboloid = Paraboloid(
      low_datum, ref_frame, aligned_paraboloid.a(), aligned_paraboloid.b());
  const auto high_paraboloid = Paraboloid(
      high_datum, ref_frame, aligned_paraboloid.a(), aligned_paraboloid.b());
  double result = getVolumeFraction(a_cell, low_paraboloid);
  if (result > target_volume_fraction_m) {
    // Volume fraction should be less than target_volume_fraction_m
    std::cout << "Lower bound result failed: " << result << std::endl;
    return false;
  } else {
    std::cout << "Lower bound result passed: " << result << std::endl;
  }
  result = getVolumeFraction(a_cell, high_paraboloid);
  if (result < target_volume_fraction_m) {
    // Volume fraction should be > target_volume_fraction_m.
    std::cout << "Upper bound result failed: " << result << std::endl;
    return false;
  } else {
    std::cout << "Upper bound result passed: " << result << std::endl;
  }
  return true;
}

template <class CellType>
void ProgressiveDistanceSolverParaboloid<CellType>::solveForDistance(
    const CellType& a_cell) {
  // Calculate volume of cell
  initial_cell_volume_m = a_cell.calculateVolume();

  // Generate half-edge structure
  auto& complete_polytope = setHalfEdgeStructureParaboloid(a_cell);
  auto half_edge_polytope =
      generateSegmentedVersionParaboloid<CellType>(&complete_polytope);
  assert(half_edge_polytope.checkValidHalfEdgeStructure());

  // Move cell to local frame of reference of the paraboloid
  const auto& datum = reconstruction_m.getDatum();
  const auto& ref_frame = reconstruction_m.getReferenceFrame();
  const auto& aligned_paraboloid = reconstruction_m.getAlignedParaboloid();
  const UnsignedIndex_t number_of_vertices =
      half_edge_polytope.getNumberOfVertices();
  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    const Pt original_pt =
        half_edge_polytope.getVertex(v)->getLocation().getPt() - datum;
    typename CellType::pt_type projected_location;
    auto& pt = projected_location.getPt();
    for (UnsignedIndex_t n = 0; n < 3; ++n) {
      pt[n] = ref_frame[n] * original_pt;
    }
    // std::cout << "New pt = " << pt << std::endl;
    half_edge_polytope.getVertex(v)->setLocation(projected_location);
  }

  // Sort the vertices Z positions in the local frame of reference
  std::vector<double> vertices_z;
  vertices_z.resize(number_of_vertices);
  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    vertices_z[v] = half_edge_polytope.getVertex(v)->getLocation().getPt()[2];
  }

  // Compute new face planes
  for (auto& face : half_edge_polytope) {
    auto normal = Normal(0.0, 0.0, 0.0);
    const auto starting_half_edge = face->getStartingHalfEdge();
    auto current_half_edge = starting_half_edge;
    auto next_half_edge = starting_half_edge->getNextHalfEdge();
    const auto& start_location =
        starting_half_edge->getPreviousVertex()->getLocation().getPt();
    do {
      normal += crossProduct(
          current_half_edge->getVertex()->getLocation().getPt() -
              start_location,
          next_half_edge->getVertex()->getLocation().getPt() - start_location);
      current_half_edge = next_half_edge;
      next_half_edge = next_half_edge->getNextHalfEdge();
    } while (next_half_edge != starting_half_edge);
    normal.normalize();
    face->setPlane(Plane(normal, normal * start_location));
  }
  assert(half_edge_polytope.checkValidHalfEdgeStructure());

  // Compute and sort distances between the paraboloid and the vertices
  // ---> in the local frame of reference of the paraboloid
  sorted_distances_m.resize(number_of_vertices + 2);
  UnsignedIndex_t count = 0;
  const Pt& start_pt = half_edge_polytope.getVertex(0)->getLocation().getPt();
  sorted_distances_m[count++] = signedDistance(start_pt, aligned_paraboloid);
  double cell_zmin = start_pt[2], cell_zmax = start_pt[2];
  for (UnsignedIndex_t v = 1; v < number_of_vertices; ++v) {
    const Pt& vertex = half_edge_polytope.getVertex(v)->getLocation().getPt();
    sorted_distances_m[count++] = signedDistance(vertex, aligned_paraboloid);
    cell_zmin = std::min({vertex[2], cell_zmin});
    cell_zmax = std::max({vertex[2], cell_zmax});
  }
  // Sort  into ascending
  sorted_distances_m[count++] = cell_zmin;
  sorted_distances_m[count++] = cell_zmax;
  std::sort(sorted_distances_m.begin(), sorted_distances_m.end());

  // std::cout << "List of sorted distances = ";
  // for (UnsignedIndex_t v = 0; v < number_of_vertices + 2; ++v) {
  //   std::cout << sorted_distances_m[v] << "  ";
  // }
  // std::cout << std::endl;
  // assert(this->isBoundsTrueBounds(a_cell));

  // Keep bisectioning between nodes until we find between where the solution
  // lays. Will then just have a prismatoid that needs to be optimized over
  std::array<UnsignedIndex_t, 3> bounding_indices{
      {0, 0, static_cast<UnsignedIndex_t>(sorted_distances_m.size()) - 1}};
  std::array<double, 3> bounding_values{{0.0, 0.0, 1.0}};
  while (bounding_indices[2] - bounding_indices[0] > 1) {
    bounding_indices[1] = (bounding_indices[0] + bounding_indices[2]) / 2;
    // std::cout << "Distance = " << sorted_distances_m[bounding_indices[1]]
    // << std::endl;
    shiftSegmentedPolytopeAndUpdateFaces(
        &half_edge_polytope, vertices_z,
        sorted_distances_m[bounding_indices[1]]);
    assert(half_edge_polytope.checkValidHalfEdgeStructure());
    // std::cout << "Caculating volume... " << std::endl;
    const double volume_fraction_cut =
        intersectPolyhedronWithParaboloid<Volume>(
            &half_edge_polytope, &complete_polytope, aligned_paraboloid) /
        initial_cell_volume_m;
    resetPolyhedron(&half_edge_polytope, &complete_polytope);
    // std::cout << "Volume calculated is " << volume_fraction_cut << std::endl;
    if (volume_fraction_cut >
        target_volume_fraction_m + volume_fraction_tolerance_m) {
      bounding_indices[2] = bounding_indices[1];
      bounding_values[2] = volume_fraction_cut;
    } else if (volume_fraction_cut <
               target_volume_fraction_m - volume_fraction_tolerance_m) {
      bounding_indices[0] = bounding_indices[1];
      bounding_values[0] = volume_fraction_cut;
    } else {
      distances_m = sorted_distances_m[bounding_indices[1]];
      return;
    }
  }

  double interval_min = sorted_distances_m[bounding_indices[0]];
  double interval_max = sorted_distances_m[bounding_indices[2]];

  // std::cout << std::endl
  //           << "Searching in the interval : ["
  //           << sorted_distances_m[bounding_indices[0]] << ", "
  //           << sorted_distances_m[bounding_indices[2]] << "]" << std::endl
  //           << std::endl;

  // // Now just have a prismatoid left. Perform secant, fall back to
  // Bisection if it fails.
  double distance = sorted_distances_m[bounding_indices[0]];
  double delta = sorted_distances_m[bounding_indices[0]] -
                 sorted_distances_m[bounding_indices[2]];
  double old_error = bounding_values[2] - target_volume_fraction_m;
  double error = bounding_values[0] - target_volume_fraction_m;
  bounding_values[0] = sorted_distances_m[bounding_indices[0]];
  bounding_values[2] = sorted_distances_m[bounding_indices[2]];
  for (UnsignedIndex_t iter = 0; iter < max_iter_m; ++iter) {
    delta *= -error / safelyEpsilon(error - old_error);
    old_error = error;
    distance += delta;
    // std::cout << "Distance = " << distance << std::endl;
    shiftSegmentedPolytopeAndUpdateFaces(&half_edge_polytope, vertices_z,
                                         distance);
    // std::cout << "Caculating volume... " << std::endl;
    const double volume_fraction_cut =
        intersectPolyhedronWithParaboloid<Volume>(
            &half_edge_polytope, &complete_polytope, aligned_paraboloid) /
        initial_cell_volume_m;
    resetPolyhedron(&half_edge_polytope, &complete_polytope);
    // std::cout << "Volume calculated is " << volume_fraction_cut << std::endl;
    error = volume_fraction_cut - target_volume_fraction_m;
    if (volume_fraction_cut >
        target_volume_fraction_m + volume_fraction_tolerance_m) {
      if (distance < bounding_values[2]) {
        bounding_values[2] = distance;
      }
    } else if (volume_fraction_cut <
               target_volume_fraction_m - volume_fraction_tolerance_m) {
      if (distance > bounding_values[0]) {
        bounding_values[0] = distance;
      }
    } else {
      distances_m = distance;
      return;
    }
  }
  bounding_values[0] = interval_min;  // std::max(bounding_values[0],
                                      // sorted_distances_m.front());
  bounding_values[2] =
      interval_max;  // std::min(bounding_values[2], sorted_distances_m.back());

  // std::cout << std::endl
  //           << "Performing bisection (fall-back)in the interval : ["
  //           << bounding_values[0] << ", " << bounding_values[2] << "]"
  //           << std::endl
  //           << std::endl;

  // Perform bisection since secant failed to find answer within tolerance.
  // Move cell back to its initial position

  for (UnsignedIndex_t iter = 0; iter < max_bisection_iter; ++iter) {
    bounding_values[1] = 0.5 * (bounding_values[0] + bounding_values[2]);
    // std::cout << "Distance =  " << bounding_values[1] << std::endl;
    shiftSegmentedPolytopeAndUpdateFaces(&half_edge_polytope, vertices_z,
                                         bounding_values[1]);
    // std::cout << "Caculating volume... " << std::endl;
    const double volume_fraction_cut =
        intersectPolyhedronWithParaboloid<Volume>(
            &half_edge_polytope, &complete_polytope, aligned_paraboloid) /
        initial_cell_volume_m;
    resetPolyhedron(&half_edge_polytope, &complete_polytope);
    // std::cout << "Volume calculated is " << volume_fraction_cut << std::endl;
    if (volume_fraction_cut >
        target_volume_fraction_m + volume_fraction_tolerance_m) {
      bounding_values[2] = bounding_values[1];
    } else if (volume_fraction_cut <
               target_volume_fraction_m - volume_fraction_tolerance_m) {
      bounding_values[0] = bounding_values[1];
    } else {
      distances_m = bounding_values[1];
      return;
    }
  }

  std::cout << "BISECTION HAS FAILED" << std::endl;
  distances_m = bounding_values[1];
  return;
}  // namespace IRL

template <class SegmentedHalfEdgePolytopeType>
void shiftSegmentedPolytopeAndUpdateFaces(
    SegmentedHalfEdgePolytopeType* a_seg_half_edge,
    const std::vector<double>& a_original_vert_Z, const double a_distance) {
  // Shift along Z
  const UnsignedIndex_t number_of_vertices =
      a_seg_half_edge->getNumberOfVertices();
  assert(a_original_vert_Z.size() == number_of_vertices);
  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    const Pt pt = a_seg_half_edge->getVertex(v)->getLocation().getPt();
    typename SegmentedHalfEdgePolytopeType::pt_type new_location;
    Pt& new_pt = new_location.getPt();
    new_pt[0] = pt[0];
    new_pt[1] = pt[1];
    new_pt[2] = a_original_vert_Z[v] - a_distance;
    // std::cout << "New point = " << new_pt << std::endl;
    a_seg_half_edge->getVertex(v)->setLocation(new_location);
  }
  // Update faces
  for (auto& face : (*a_seg_half_edge)) {
    const auto normal = face->getPlane().normal();
    const auto starting_half_edge = face->getStartingHalfEdge();
    const auto& start_location =
        starting_half_edge->getPreviousVertex()->getLocation().getPt();
    // std::cout << "Pt ref = " << start_location << std::endl;
    // std::cout << "New distance = " << normal * start_location << std::endl;
    face->setPlane(Plane(normal, normal * start_location));
  }
}

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PROGRESSIVE_DISTANCE_SOLVER_PARABOLOID_TPP_
