// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYGON_TPP_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYGON_TPP_

namespace IRL {

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
SegmentedHalfEdgePolygon<FaceType, VertexType>
HalfEdgePolygon<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                kMaxVertices, kMaxFaces>::generateSegmentedPolygon(void) {
  SegmentedHalfEdgePolygon<FaceType, VertexType> a_polygon;
  this->setSegmentedPolygon(&a_polygon);
  return a_polygon;
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
template <class SegmentedType>
void HalfEdgePolygon<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::setSegmentedPolygon(SegmentedType* a_polytope) {
  a_polytope->setNumberOfFaces(0);
  a_polytope->setNumberOfVertices(0);
  for (UnsignedIndex_t n = 0; n < this->getNumberOfFaces(); ++n) {
    a_polytope->addFace(&this->getFace(n));
  }
  for (UnsignedIndex_t n = 0; n < this->getNumberOfVertices(); ++n) {
    a_polytope->addVertex(&this->getVertex(n));
  }
  a_polytope->setPlaneOfExistence(&this->getPlaneOfExistence());
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
void HalfEdgePolygon<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                     kMaxVertices,
                     kMaxFaces>::setPlaneOfExistence(const Plane& a_plane) {
  plane_of_existence_m = a_plane;
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
const Plane&
HalfEdgePolygon<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                kMaxVertices, kMaxFaces>::getPlaneOfExistence(void) {
  return plane_of_existence_m;
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYGON_TPP_
