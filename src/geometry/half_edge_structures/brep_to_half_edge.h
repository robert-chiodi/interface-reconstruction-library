// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_BREP_TO_HALF_EDGE_H_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_BREP_TO_HALF_EDGE_H_

#include "src/parameters/defined_types.h"
#include "src/geometry/half_edge_structures/half_edge_polyhedron.h"
#include "src/geometry/half_edge_structures/half_edge.h"

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
class BREPToHalfEdge : public HalfEdgePolytope<PtType> {
 public:

  BREPToHalfEdge(void) = delete;

  template <class PtContainerType, class FaceBREPType>
    HalfEdgePolyhedron<PtType, VertexType,HalfEdgeType,FaceType,kMaxHalfEdges, kMaxVertices, kMaxFaces>
  static generateHalfEdgeVersion(const PtContainerType& a_pt_container, const FaceBREPType& a_face_brep);

  template <class PtContainerType, class FaceBREPType, class HalfEdgePolyhedronType>
  static void setHalfEdgeVersion(const PtContainerType& a_pt_container,
			    const FaceBREPType& a_face_brep,
			    HalfEdgePolyhedronType* a_half_edge_version);
  
  ~BREPToHalfEdge(void) = delete;

 private:
};
}  // namespace IRL

#include "src/geometry/half_edge_structures/brep_to_half_edge.tpp"

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_BREP_TO_HALF_EDGE_H_
