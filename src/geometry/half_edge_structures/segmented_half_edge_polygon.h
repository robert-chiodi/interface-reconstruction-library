// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYGON_H_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYGON_H_

#include <unordered_map>

#include "src/geometry/half_edge_structures/half_edge.h"
#include "src/geometry/half_edge_structures/segmented_half_edge_polytope.h"
#include "src/moments/volume.h"
#include "src/moments/volume_moments.h"
#include "src/moments/volume_moments_and_normal.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <
    class FaceType, class VertexType,
    UnsignedIndex_t kMaxFaces =
        segmented_half_edge_polytope::default_sizes::segemented_max_faces,
    UnsignedIndex_t kMaxVertices =
        segmented_half_edge_polytope::default_sizes::segmented_max_vertices>
class SegmentedHalfEdgePolygonCommon
    : public SegmentedHalfEdgePolytope<FaceType, VertexType> {
 public:
  // This plane of existence should be the plane that the
  // polygon lives on. It is used to assign a consistent
  // reference around which to build signed moments during
  // calculateMoments. Most likely, this should be set during
  // generation of the SegmentedHalfEdgePolygonCommon, perhaps from
  // a HalfEdgePolygon generated from another Polygon class,
  // or during a subsequent splitting operation, where it will
  // assume its plane of existence from the SegmentedHalfEdgePolygonCommon
  // it was split from (which presumably, going back far enough, was initialized
  // during the original creation of the first, over-arching, HalfEdgePolygon
  // from a different Polygon class).
  void setPlaneOfExistence(const Plane* a_plane);
  const Plane& getPlaneOfExistence(void) const;

  inline Volume calculateVolume(void) const;

  inline Volume calculateAbsoluteVolume(void) const;

  inline VolumeMomentsAndNormal calculateVolumeMomentsAndNormal(void) const;

  inline Pt calculateCentroid(void) const;

  inline VolumeMoments calculateMoments(void) const;

  bool checkValidHalfEdgeStructure(void);

 private:
  bool checkVerticesCompleteCycle(void) const;
  bool checkVerticesOnlyHaveHalfEdgesWithCorrectFaces(void) const;

  const Plane* plane_of_existence_m;
};

template <
    class PtType, class FaceType, class VertexType,
    UnsignedIndex_t kMaxFaces =
        segmented_half_edge_polytope::default_sizes::segemented_max_faces,
    UnsignedIndex_t kMaxVertices =
        segmented_half_edge_polytope::default_sizes::segmented_max_vertices>
class SegmentedHalfEdgePolygonSpecificPt;

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
class SegmentedHalfEdgePolygonSpecificPt<Pt, FaceType, VertexType, kMaxFaces,
                                         kMaxVertices>
    : public SegmentedHalfEdgePolygonCommon<FaceType, VertexType, kMaxFaces,
                                            kMaxVertices> {
 public:
  using SegmentedHalfEdgePolygonCommon<
      FaceType, VertexType, kMaxFaces,
      kMaxVertices>::SegmentedHalfEdgePolygonCommon;
};

template <
    class FaceType, class VertexType,
    UnsignedIndex_t kMaxFaces =
        segmented_half_edge_polytope::default_sizes::segemented_max_faces,
    UnsignedIndex_t kMaxVertices =
        segmented_half_edge_polytope::default_sizes::segmented_max_vertices>
class SegmentedHalfEdgePolygon
    : public SegmentedHalfEdgePolygonSpecificPt<typename VertexType::pt_type,
                                                FaceType, VertexType, kMaxFaces,
                                                kMaxVertices> {
 public:
  using SegmentedHalfEdgePolygonSpecificPt<
      typename VertexType::pt_type, FaceType, VertexType, kMaxFaces,
      kMaxVertices>::SegmentedHalfEdgePolygonSpecificPt;
};

}  // namespace IRL

#include "src/geometry/half_edge_structures/segmented_half_edge_polygon.tpp"

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYGON_H_
