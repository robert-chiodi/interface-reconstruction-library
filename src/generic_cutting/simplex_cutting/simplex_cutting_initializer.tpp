// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_SIMPLEX_CUTTING_SIMPLEX_CUTTING_INITIALIZER_TPP_
#define SRC_GENERIC_CUTTING_SIMPLEX_CUTTING_SIMPLEX_CUTTING_INITIALIZER_TPP_

#include "src/geometry/decomposed_polytope/decomposed_polytope_vertex_storage.h"
#include "src/geometry/decomposed_polytope/segmented_decomposed_polytope.h"
#include "src/geometry/general/geometry_type_traits.h"
#include "src/geometry/polyhedrons/tet.h"
#include "src/helpers/SFINAE_boiler_plate.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class EncompassingGeometryType>
enable_if_t<
    is_polygon<EncompassingGeometryType>::value,
    DecomposedPolygonVertexStorage<typename EncompassingGeometryType::pt_type>&>
getVertexStorage(const EncompassingGeometryType& a_geometry);

template <class EncompassingGeometryType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value,
            DecomposedPolyhedronVertexStorage<
                typename EncompassingGeometryType::pt_type>&>
getVertexStorage(const EncompassingGeometryType& a_geometry);

template <class EncompassingGeometryType, class VertexStorageType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value,
            SegmentedDecomposedPolyhedron<typename VertexStorageType::pt_type>>
getInitialSimplexList(const EncompassingGeometryType& a_geometry,
                      VertexStorageType* a_vertex_storage);

template <class EncompassingGeometryType, class VertexStorageType>
enable_if_t<is_polygon<EncompassingGeometryType>::value,
            SegmentedDecomposedPolygon<typename VertexStorageType::pt_type>>
getInitialSimplexList(const EncompassingGeometryType& a_geometry,
                      VertexStorageType* a_vertex_storage);

//******************************************************************* //
//     Template function definitions placed below this.
//******************************************************************* //
template <class ReturnType, class EncompassingType, class ReconstructionType>
inline ReturnType cutThroughSimplex(
    const EncompassingType& a_encompassing_polyhedron,
    const ReconstructionType& a_separating_reconstruction) {
  auto& vertex_storage = getVertexStorage(a_encompassing_polyhedron);
  auto simplex_list =
      getInitialSimplexList(a_encompassing_polyhedron, &vertex_storage);
  return getVolumeMoments<ReturnType, SimplexCutting>(
      &simplex_list, &vertex_storage, a_separating_reconstruction);
}

// Static storage for each kind of vertex type in polygon.
template <class EncompassingGeometryType>
enable_if_t<
    is_polygon<EncompassingGeometryType>::value,
    DecomposedPolygonVertexStorage<typename EncompassingGeometryType::pt_type>&>
getVertexStorage(const EncompassingGeometryType& a_geometry) {
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
                typename EncompassingGeometryType::pt_type>&>
getVertexStorage(const EncompassingGeometryType& a_geometry) {
  static DecomposedPolyhedronVertexStorage<
      typename EncompassingGeometryType::pt_type>
      vertex_storage;
  vertex_storage.resetFromGeometry(a_geometry);
  return vertex_storage;
}

template <class EncompassingGeometryType, class VertexStorageType>
enable_if_t<is_polyhedron<EncompassingGeometryType>::value,
            SegmentedDecomposedPolyhedron<typename VertexStorageType::pt_type>>
getInitialSimplexList(const EncompassingGeometryType& a_geometry,
                      VertexStorageType* a_vertex_storage) {
  return SegmentedDecomposedPolyhedron<typename VertexStorageType::pt_type>(
      a_geometry, *a_vertex_storage);
}

template <class EncompassingGeometryType, class VertexStorageType>
enable_if_t<is_polygon<EncompassingGeometryType>::value,
            SegmentedDecomposedPolygon<typename VertexStorageType::pt_type>>
getInitialSimplexList(const EncompassingGeometryType& a_geometry,
                      VertexStorageType* a_vertex_storage) {
  SegmentedDecomposedPolygon<typename VertexStorageType::pt_type>
      segmented_polygon(a_geometry, *a_vertex_storage);
  segmented_polygon.setPlaneOfExistence(
      &a_vertex_storage->getPlaneOfExistence());
  return segmented_polygon;
}

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_SIMPLEX_CUTTING_SIMPLEX_CUTTING_INITIALIZER_TPP_
