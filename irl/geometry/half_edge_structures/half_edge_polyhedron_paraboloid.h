// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYHEDRON_PARABOLOID_H_
#define IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYHEDRON_PARABOLOID_H_

#include "irl/geometry/half_edge_structures/half_edge_paraboloid.h"
#include "irl/geometry/half_edge_structures/half_edge_polytope.h"
#include "irl/geometry/half_edge_structures/segmented_half_edge_polyhedron_paraboloid.h"
#include "irl/parameters/defined_types.h"

namespace IRL {
template <class PtType, class VertexType = VertexParaboloid<PtType>,
          class HalfEdgeType = HalfEdgeParaboloid<VertexType>,
          class FaceType = FaceParaboloid<HalfEdgeType>,
          UnsignedIndex_t kMaxHalfEdges =
              half_edge_polytope::default_sizes::complete_block_half_edges,
          UnsignedIndex_t kMaxVertices =
              half_edge_polytope::default_sizes::complete_block_vertices,
          UnsignedIndex_t kMaxFaces =
              half_edge_polytope::default_sizes::complete_block_faces>
class HalfEdgePolyhedronParaboloid
    : public HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType,
                              kMaxHalfEdges, kMaxVertices, kMaxFaces> {
 public:
  HalfEdgePolyhedronParaboloid(void) = default;

  SegmentedHalfEdgePolyhedronParaboloid<FaceType, VertexType>
  generateSegmentedPolyhedron(void);

  template <class SegmentedType>
  void setSegmentedPolyhedron(SegmentedType* a_polytope);

  ~HalfEdgePolyhedronParaboloid(void) = default;

 private:
};
}  // namespace IRL

#include "irl/geometry/half_edge_structures/half_edge_polyhedron_paraboloid.tpp"

#endif  // IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYHEDRON_PARABOLOID_H_
