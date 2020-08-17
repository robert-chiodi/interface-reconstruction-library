// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYHEDRON_H_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYHEDRON_H_

#include "src/geometry/half_edge_structures/half_edge_polytope.h"
#include "src/geometry/half_edge_structures/segmented_half_edge_polyhedron.h"
#include "src/parameters/defined_types.h"

namespace IRL {
template <class PtType, class VertexType = Vertex<PtType>,
          class HalfEdgeType = HalfEdge<VertexType>,
          class FaceType = Face<HalfEdgeType>,
          UnsignedIndex_t kMaxHalfEdges =
              half_edge_polytope::default_sizes::complete_block_half_edges,
          UnsignedIndex_t kMaxVertices =
              half_edge_polytope::default_sizes::complete_block_vertices,
          UnsignedIndex_t kMaxFaces =
              half_edge_polytope::default_sizes::complete_block_faces>
class HalfEdgePolyhedron : public HalfEdgePolytope<PtType> {
 public:
  HalfEdgePolyhedron(void) = default;

  SegmentedHalfEdgePolyhedron<FaceType, VertexType> generateSegmentedPolyhedron(
      void);

  template <class SegmentedType>
  void setSegmentedPolyhedron(SegmentedType* a_polytope);

  ~HalfEdgePolyhedron(void) = default;

 private:
};
}  // namespace IRL

#include "src/geometry/half_edge_structures/half_edge_polyhedron.tpp"

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYHEDRON_H_
