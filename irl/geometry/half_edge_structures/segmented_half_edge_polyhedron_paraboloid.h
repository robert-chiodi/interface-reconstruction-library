// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYHEDRON_PARABOLOID_H_
#define IRL_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYHEDRON_PARABOLOID_H_

#include <utility>

#include "irl/data_structures/small_vector.h"
#include "irl/geometry/general/geometry_type_traits.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/half_edge_structures/half_edge.h"
#include "irl/geometry/half_edge_structures/segmented_half_edge_polyhedron.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

namespace segmented_half_edge_polyhedron_paraboloid_detail {
template <class GeometryType, class CalculationFunctor>
auto calculateMoments(GeometryType* a_geometry,
                      CalculationFunctor a_moment_accumulator) ->
    typename CalculationFunctor::ReturnType;

template <class VertexType>
class TetPtReferenceWrapper {
 public:
  using pt_type = VertexType;

  TetPtReferenceWrapper(void) = default;

  const VertexType& operator[](const UnsignedIndex_t a_index) const {
    assert(a_index < 4);
    assert(tet_m[a_index] != nullptr);
    return *tet_m[a_index];
  }

  void setPt(const UnsignedIndex_t a_index, const VertexType* a_pt) {
    assert(a_index < 4);
    assert(a_pt != nullptr);
    tet_m[a_index] = a_pt;
  }

  ~TetPtReferenceWrapper(void) = default;

 private:
  std::array<const VertexType*, 4> tet_m;
};

template <class GeometryType, class CalculationFunctor>
auto calculateMoments(GeometryType* a_geometry,
                      CalculationFunctor a_moment_accumulator) ->
    typename CalculationFunctor::ReturnType {
  TetPtReferenceWrapper<typename GeometryType::pt_type> facet_and_pt;
  // Mark all faces that should be included as not visited
  for (auto& face : *a_geometry) {
    face->markAsNotVisited();
  }

  // Now use common point (arbitrarily) as the first vertex. Cycle around and
  // mark all faces visited.
  auto& datum_vertex = *(a_geometry->getVertex(0));
  facet_and_pt.setPt(3, &datum_vertex.getLocation());
  auto current_half_edge = datum_vertex.getHalfEdge();
  do {
    current_half_edge->getFace()->markAsVisited();
    current_half_edge =
        current_half_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
  } while (current_half_edge != datum_vertex.getHalfEdge());

  // Now loop over faces that haven't been visited, triangulate, form tets,
  // and sum up moments.
  for (auto& face : (*a_geometry)) {
    if (face->hasNotBeenVisited()) {
      facet_and_pt.setPt(
          0, &face->getStartingHalfEdge()->getVertex()->getLocation());
      auto first_half_edge = face->getStartingHalfEdge()->getNextHalfEdge();
      auto second_half_edge = first_half_edge->getNextHalfEdge();
      do {
        facet_and_pt.setPt(1, &first_half_edge->getVertex()->getLocation());
        facet_and_pt.setPt(2, &second_half_edge->getVertex()->getLocation());
        a_moment_accumulator(facet_and_pt);

        first_half_edge = second_half_edge;
        second_half_edge = second_half_edge->getNextHalfEdge();
      } while (second_half_edge != face->getStartingHalfEdge());
    }
  }
  return a_moment_accumulator.getMoments();
}

}  // namespace segmented_half_edge_polyhedron_paraboloid_detail

template <
    class FaceType, class VertexType,
    UnsignedIndex_t kMaxFaces =
        segmented_half_edge_polytope::default_sizes::segemented_max_faces,
    UnsignedIndex_t kMaxVertices =
        segmented_half_edge_polytope::default_sizes::segmented_max_vertices>
class SegmentedHalfEdgePolyhedronParaboloid
    : public SegmentedHalfEdgePolyhedron<FaceType, VertexType, kMaxFaces,
                                         kMaxVertices> {
  using Base = SegmentedHalfEdgePolyhedron<FaceType, VertexType, kMaxFaces,
                                           kMaxVertices>;

 public:
  Volume calculateVolume(void) {
    return segmented_half_edge_polyhedron_paraboloid_detail::calculateMoments(
        this, Volume3D_Functor());
  }

 private:
};

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
struct is_polyhedron<SegmentedHalfEdgePolyhedronParaboloid<
    FaceType, VertexType, kMaxFaces, kMaxVertices>> : std::true_type {};

}  // namespace IRL

#endif  // IRL_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYHEDRON_PARABOLOID_H_
