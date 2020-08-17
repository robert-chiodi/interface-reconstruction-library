// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_PLANE_DISTANCE_TPP_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_PLANE_DISTANCE_TPP_

#include "src/generic_cutting/simplex_cutting/simplex_cutting_initializer.h"
#include "src/geometry/decomposed_polytope/decomposed_polytope_vertex_storage.h"
#include "src/geometry/decomposed_polytope/segmented_decomposed_polytope.h"
#include "src/helpers/SFINAE_boiler_plate.h"
#include "src/interface_reconstruction_methods/progressive_distance_solver.h"

namespace IRL {

template <class CellType, class ReconstructionType>
inline void
runIterativeSolverForDistance(const CellType &a_cell,
                              const double a_volume_fraction,
                              ReconstructionType *a_reconstruction,
                              const double a_volume_fraction_tolerance) {
  assert(a_reconstruction != nullptr);
  IterativeSolverForDistance<ReconstructionDefaultCuttingMethod, CellType>
      solver(a_cell, a_volume_fraction, a_volume_fraction_tolerance,
             (*a_reconstruction));
  a_reconstruction->setDistances(solver.getDistances());
}

template <class CellType, class ReconstructionType>
inline void
runProgressiveDistanceSolver(const CellType &a_cell,
                             const double a_volume_fraction,
                             ReconstructionType *a_reconstruction,
                             const double a_volume_fraction_tolerance) {
  assert(a_reconstruction != nullptr);
  assert(a_reconstruction->getNumberOfPlanes() == 1);
  ProgressiveDistanceSolver<CellType> solver(a_cell, a_volume_fraction,
                                             a_volume_fraction_tolerance,
                                             (*a_reconstruction));
  (*a_reconstruction)[0].distance() = solver.getDistances(0);
}

namespace plane_distance_details {
template <class VolumeFractionArrayType>
static inline double
sumVolumeFraction(const VolumeFractionArrayType &a_volume_fraction) {
  double sum = 0.0;
  for (const auto &elem : a_volume_fraction) {
    sum += elem;
  }
  return sum;
}

template <class EncompassingGeometryType>
enable_if_t<is_polygon<EncompassingGeometryType>::value,
            DecomposedPolygonVertexStorage<
                typename EncompassingGeometryType::pt_type> &>
getVertexStorage(const EncompassingGeometryType &a_geometry) {
  static DecomposedPolygonVertexStorage<
      typename EncompassingGeometryType::pt_type>
      vertex_storage;
  vertex_storage.resetFromGeometry(a_geometry);
  vertex_storage.setPlaneOfExistence(a_geometry.getPlaneOfExistence());
  return vertex_storage;
}

// Static storage for each kind of vertex type in polyhedron.
template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value,
            DecomposedPolyhedronVertexStorage<
                typename EncompassingGeometryType::pt_type> &>
getVertexStorage(const EncompassingGeometryType &a_geometry) {
  static DecomposedPolyhedronVertexStorage<
      typename EncompassingGeometryType::pt_type>
      vertex_storage;
  vertex_storage.resetFromGeometry(a_geometry);
  return vertex_storage;
}
} // namespace plane_distance_details

template <class CellType, class VolumeFractionArrayType>
inline void
runIterativeSolverForDistance(const CellType &a_cell,
                              const VolumeFractionArrayType &a_volume_fraction,
                              PlanarSeparatorPathGroup *a_reconstruction,
                              const double a_volume_fraction_tolerance) {
  assert(a_reconstruction != nullptr);
  assert(fabs(plane_distance_details::sumVolumeFraction(a_volume_fraction) -
              1.0) < 10.0 * DBL_EPSILON);
  assert(a_reconstruction->getPriorityOrderSize() == a_volume_fraction.size());
  auto current_reconstruction = a_reconstruction->getFirstReconstruction();

  // Generate a SegmentedDecomposedPolytope
  auto &vertex_storage = plane_distance_details::getVertexStorage(a_cell);
  auto split_cell = getInitialSimplexList(
      a_cell,
      &vertex_storage); // From
                        // generic_cutting/simplex_cutting/simplex_cutting_initializer.tpp

  UnsignedIndex_t vf_index = 0;
  double previous_vf_sum = 0.0;
  while (true) {
    if (!current_reconstruction.hasNeighbor()) {
      break; // Last reconstruction in chain. Break and set as all
             // encompassing.
    }
    auto &underlying_reconstruction =
        current_reconstruction.getCurrentReconstruction();
    assert(underlying_reconstruction.getNumberOfPlanes() ==
           1); // Require single plane reconstructions for now.

    const Pt centroid = split_cell.calculateCentroid();
    underlying_reconstruction[0].distance() =
        underlying_reconstruction[0].normal() * centroid;

    // Iterate for distance solution then store.
    const double scaled_volume_fraction =
        a_volume_fraction[vf_index] / (1.0 - previous_vf_sum);
    IterativeSolverForDistance<RecursiveSimplexCutting, decltype(split_cell)>
        solver(split_cell, scaled_volume_fraction, a_volume_fraction_tolerance,
               underlying_reconstruction);
    underlying_reconstruction.setDistances(solver.getDistances());

    // Truncate with flipped plane, giving only part above plane.
    const auto flipped_plane =
        underlying_reconstruction[0].generateFlippedPlane();
    truncateDecomposedPolytope(&split_cell, &vertex_storage, flipped_plane);

    previous_vf_sum += a_volume_fraction[vf_index];
    ++vf_index;
    current_reconstruction = current_reconstruction.getNeighbor();
  }
  auto &underlying_reconstruction =
      current_reconstruction.getCurrentReconstruction();
  underlying_reconstruction.setNumberOfPlanes(1);
  underlying_reconstruction[0] = Plane(Normal(0.0, 0.0, 0.0), 1.0);
}

} // namespace IRL

#endif // SRC_INTERFACE_RECONSTRUCTION_METHODS_PLANE_DISTANCE_TPP_
