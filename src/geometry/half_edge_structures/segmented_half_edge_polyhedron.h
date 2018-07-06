// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYHEDRON_H_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYHEDRON_H_

#include <unordered_map>

#include "src/geometry/general/pt.h"
#include "src/geometry/general/pt_with_data.h"
#include "src/geometry/half_edge_structures/half_edge.h"
#include "src/geometry/half_edge_structures/segmented_half_edge_polytope.h"
#include "src/moments/volume.h"
#include "src/moments/volume_moments.h"
#include "src/moments/volume_moments_and_doubles.h"
#include "src/parameters/defined_types.h"

namespace IRL {
template <
    class FaceType, class VertexType,
    UnsignedIndex_t kMaxFaces =
        segmented_half_edge_polytope::default_sizes::segemented_max_faces,
    UnsignedIndex_t kMaxVertices =
        segmented_half_edge_polytope::default_sizes::segmented_max_vertices>
class SegmentedHalfEdgePolyhedronCommon
    : public SegmentedHalfEdgePolytope<FaceType, VertexType> {
 public:
  inline Volume calculateVolume(void);

  inline Volume calculateAbsoluteVolume(void);

  inline Pt calculateCentroid(void);

  inline VolumeMoments calculateMoments(void);

  bool checkValidHalfEdgeStructure(void);

 private:
  void markFacesTouchedByOpposite(
      FaceType* a_face_to_touch_from,
      std::unordered_map<FaceType*, bool>* a_map_to_fill);
  bool checkVerticesCompleteCycle(void) const;
  bool checkVerticesOnlyHaveHalfEdgesWithCorrectFaces(void) const;
};

template <
    class PtType, class FaceType, class VertexType,
    UnsignedIndex_t kMaxFaces =
        segmented_half_edge_polytope::default_sizes::segemented_max_faces,
    UnsignedIndex_t kMaxVertices =
        segmented_half_edge_polytope::default_sizes::segmented_max_vertices>
class SegmentedHalfEdgePolyhedronSpecificPt
    : public SegmentedHalfEdgePolyhedronCommon<FaceType, VertexType, kMaxFaces,
                                               kMaxVertices> {};

template <class FunctorType, UnsignedIndex_t kArrayLength, class FaceType,
          class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
class SegmentedHalfEdgePolyhedronSpecificPt<
    PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>, FaceType,
    VertexType, kMaxFaces, kMaxVertices>
    : public SegmentedHalfEdgePolyhedronCommon<FaceType, VertexType, kMaxFaces,
                                               kMaxVertices> {
 public:
  inline VolumeMomentsAndDoubles<kArrayLength> calculateVolumeMomentsAndDoubles(
      void);
};

template <
    class FaceType, class VertexType,
    UnsignedIndex_t kMaxFaces =
        segmented_half_edge_polytope::default_sizes::segemented_max_faces,
    UnsignedIndex_t kMaxVertices =
        segmented_half_edge_polytope::default_sizes::segmented_max_vertices>
class SegmentedHalfEdgePolyhedron
    : public SegmentedHalfEdgePolyhedronSpecificPt<typename VertexType::pt_type,
                                                   FaceType, VertexType,
                                                   kMaxFaces, kMaxVertices> {
 public:
  using SegmentedHalfEdgePolyhedronSpecificPt<
      typename VertexType::pt_type, FaceType, VertexType, kMaxFaces,
      kMaxVertices>::SegmentedHalfEdgePolyhedronSpecificPt;
};

}  // namespace IRL

#include "src/geometry/half_edge_structures/segmented_half_edge_polyhedron.tpp"

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYHEDRON_H_
