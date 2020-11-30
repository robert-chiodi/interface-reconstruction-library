// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYTOPE_H_
#define IRL_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYTOPE_H_

#include <array>

#include "irl/data_structures/small_vector.h"
#include "irl/geometry/half_edge_structures/half_edge.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

namespace segmented_half_edge_polytope {
namespace default_sizes {
static constexpr UnsignedIndex_t segemented_max_faces = 64;
static constexpr UnsignedIndex_t segmented_max_vertices = 64;
}  // namespace default_sizes
}  // namespace segmented_half_edge_polytope

// Holds segmented portions of a HalfEdgePolytope. Only includes faces! Need to
// be careful that all half-edges and vertices pointed to by these faces remain
// valid while this is being used. This allows easy linked cutting though by
// generating new SegmentedHalfEdgePolytopes during cutting, keeping an accurate
// list of the faces that make the polyhedron that represent the remaining
// polyhedron.
template <
    class FaceType, class VertexType,
    UnsignedIndex_t kMaxFaces =
        segmented_half_edge_polytope::default_sizes::segemented_max_faces,
    UnsignedIndex_t kMaxVertices =
        segmented_half_edge_polytope::default_sizes::segmented_max_vertices>
class SegmentedHalfEdgePolytope {
  using iterator = typename SmallVector<FaceType*, kMaxFaces>::iterator;
  using const_iterator =
      typename SmallVector<FaceType*, kMaxFaces>::const_iterator;

 public:
  using face_type = FaceType;
  using vertex_type = VertexType;
  using pt_type = typename VertexType::pt_type;
  static constexpr UnsignedIndex_t maxVertices = kMaxVertices;
  static constexpr UnsignedIndex_t maxFaces = kMaxFaces;

  SegmentedHalfEdgePolytope(void) = default;

  void setNumberOfFaces(const UnsignedIndex_t a_number_of_faces);
  void setNumberOfVertices(const UnsignedIndex_t a_size);

  UnsignedIndex_t getNumberOfFaces(void) const;
  UnsignedIndex_t getNumberOfVertices(void) const;

  FaceType* operator[](const UnsignedIndex_t a_index);

  const FaceType* operator[](const UnsignedIndex_t a_index) const;

  void addFace(FaceType* a_face);
  FaceType*& getFacePointer(const UnsignedIndex_t a_index);

  void addVertex(VertexType* a_vertex);
  VertexType* getVertex(const UnsignedIndex_t a_index);
  const VertexType* getVertex(const UnsignedIndex_t a_index) const;
  VertexType*& getVertexPointer(const UnsignedIndex_t a_index);

  /// \brief Return a point for the lower limits of the polyhedron in 3D space.
  IRL::Pt getLowerLimits(void) const;
  /// \brief Return a point for the upper limits of the polyhedron in 3D space.
  IRL::Pt getUpperLimits(void) const;
  /// \brief Return bounding box with lower and upper points (at index 0 and 1
  /// respectively)
  std::array<Pt, 2> getBoundingBox(void) const;

  int calculateAndStoreDistanceToVertices(const Plane& a_plane);
  void removeFace(const UnsignedIndex_t a_index);
  void removeVertex(const UnsignedIndex_t a_index);

  template <class HalfEdgePolytopeType>
  void clear(HalfEdgePolytopeType* a_complete_polytope);

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  ~SegmentedHalfEdgePolytope(void) = default;

 private:
  void checkIfStaticAllocationExceeded(void) const;

  SmallVector<FaceType*, kMaxFaces> faces_m;
  SmallVector<VertexType*, kMaxVertices> vertices_m;
};

// IO
template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
inline std::ostream& operator<<(
    std::ostream& out,
    const SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                                    kMaxVertices>& a_polyhedron);

}  // namespace IRL

#include "irl/geometry/half_edge_structures/segmented_half_edge_polytope.tpp"

#endif // IRL_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYTOPE_H_
