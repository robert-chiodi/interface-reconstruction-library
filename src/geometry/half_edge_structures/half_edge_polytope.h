// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYTOPE_H_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYTOPE_H_

#include <utility>

#include "src/data_structures/chained_block_storage.h"
#include "src/geometry/general/pt.h"
#include "src/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace half_edge_polytope {
namespace default_sizes {
static constexpr UnsignedIndex_t complete_block_vertices = 64;
static constexpr UnsignedIndex_t complete_block_faces = 32;
static constexpr UnsignedIndex_t complete_block_half_edges = 128;
}  // namespace default_sizes
}  // namespace half_edge_polytope

template <class PtType, class VertexType = Vertex<PtType>,
          class HalfEdgeType = HalfEdge<VertexType>,
          class FaceType = Face<HalfEdgeType>,
          UnsignedIndex_t kMaxHalfEdges =
              half_edge_polytope::default_sizes::complete_block_half_edges,
          UnsignedIndex_t kMaxVertices =
              half_edge_polytope::default_sizes::complete_block_vertices,
          UnsignedIndex_t kMaxFaces =
              half_edge_polytope::default_sizes::complete_block_faces>
class HalfEdgePolytope {
 public:
  using pt_type = PtType;
  using vertex_type = VertexType;
  using half_edge_type = HalfEdgeType;
  using face_type = FaceType;

  static constexpr UnsignedIndex_t maxHalfEdges = kMaxHalfEdges;
  static constexpr UnsignedIndex_t maxVertices = kMaxVertices;
  static constexpr UnsignedIndex_t maxFaces = kMaxFaces;

  HalfEdgePolytope(void) = default;

  static HalfEdgePolytope fromKnownSizes(
      const UnsignedIndex_t a_number_of_half_edges,
      const UnsignedIndex_t a_number_of_vertices,
      const UnsignedIndex_t a_number_of_faces);

  void reset(void);

  void resize(const UnsignedIndex_t a_number_of_half_edges,
              const UnsignedIndex_t a_number_of_vertices,
              const UnsignedIndex_t a_number_of_faces);

  UnsignedIndex_t getNumberOfFaces(void) const;
  UnsignedIndex_t getNumberOfVertices(void) const;

  HalfEdgeType& getHalfEdge(const UnsignedIndex_t a_index);

  const HalfEdgeType& getHalfEdge(const UnsignedIndex_t a_index) const;

  VertexType& getVertex(const UnsignedIndex_t a_index);
  const VertexType& getVertex(const UnsignedIndex_t a_index) const;
  FaceType& getFace(const UnsignedIndex_t a_index);
  const FaceType& getFace(const UnsignedIndex_t a_index) const;

  HalfEdgeType* getNewHalfEdge(void);
  HalfEdgeType* getNewHalfEdge(const HalfEdgeType& a_half_edge);
  HalfEdgeType* getNewHalfEdge(HalfEdgeType&& a_half_edge);

  VertexType* getNewVertex(void);
  VertexType* getNewVertex(VertexType&& a_vertex);

  FaceType* getNewFace(void);
  FaceType* getNewFace(FaceType&& a_face);

  template <class GeometryType>
  void setVertexLocations(const GeometryType& a_geometry);

 private:
  HalfEdgePolytope(const UnsignedIndex_t a_number_of_half_edges,
                   const UnsignedIndex_t a_number_of_vertices,
                   const UnsignedIndex_t a_number_of_faces);

  ChainedBlockStorage<HalfEdgeType, kMaxHalfEdges> half_edges_m;
  ChainedBlockStorage<VertexType, kMaxVertices> vertices_m;
  ChainedBlockStorage<FaceType, kMaxFaces> faces_m;
};

// IO
template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
inline std::ostream& operator<<(
    std::ostream& out,
    const HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType,
                           kMaxHalfEdges, kMaxVertices, kMaxFaces>&
        a_polyhedron);
}  // namespace IRL

#include "src/geometry/half_edge_structures/half_edge_polytope.tpp"

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYTOPE_H_
