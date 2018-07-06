// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYGON_TPP_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYGON_TPP_

#include "src/geometry/general/moment_calculation_through_simplices.h"

namespace IRL {

namespace segmented_half_edge_polygon_detail {
template <class GeometryType, class CalculationFunctor>
auto calculateMoments(GeometryType* a_geometry,
                      CalculationFunctor a_moment_accumulator) ->
    typename CalculationFunctor::ReturnType;

template <class VertexType>
class TriPtReferenceWrapper {
 public:
  using pt_type = VertexType;

  TriPtReferenceWrapper(void) = delete;

  TriPtReferenceWrapper(const Plane& a_plane_of_existence)
      : plane_of_existence_m(a_plane_of_existence) {}

  const VertexType& operator[](const UnsignedIndex_t a_index) const {
    assert(a_index < 3);
    assert(triangle_m[a_index] != nullptr);
    return *triangle_m[a_index];
  }

  void setPt(const UnsignedIndex_t a_index, const VertexType* a_pt) {
    assert(a_index < 3);
    assert(a_pt != nullptr);
    triangle_m[a_index] = a_pt;
  }

  const Plane& getPlaneOfExistence(void) const { return plane_of_existence_m; }

  ~TriPtReferenceWrapper(void) = default;

 private:
  const Plane& plane_of_existence_m;
  std::array<const VertexType*, 3> triangle_m;
};

template <class GeometryType, class CalculationFunctor>
auto calculateMoments(GeometryType* a_geometry,
                      CalculationFunctor a_moment_accumulator) ->
    typename CalculationFunctor::ReturnType {

  if(a_geometry->getNumberOfFaces() == 0){
	return a_moment_accumulator.getMoments();
  }

  segmented_half_edge_polygon_detail::TriPtReferenceWrapper<
      typename GeometryType::pt_type>
      facet(a_geometry->getPlaneOfExistence());
  facet.setPt(0, &(a_geometry->getVertex(0)->getLocation()));

  // Now loop over faces that haven't been visited, triangulate, form tets,
  // and sum up moments.
  for (auto& face : (*a_geometry)) {
    if (face == &getOpenBoundaryFace<typename GeometryType::face_type>()) {
      continue;
    }
    auto current_half_edge = face->getStartingHalfEdge();
    do {
      facet.setPt(1, &(current_half_edge->getPreviousVertex()->getLocation()));
      facet.setPt(2, &(current_half_edge->getVertex()->getLocation()));
      a_moment_accumulator(facet);
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != face->getStartingHalfEdge());
  }
  return a_moment_accumulator.getMoments();
}
}  // namespace segmented_half_edge_polygon_detail

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
void SegmentedHalfEdgePolygonCommon<
    FaceType, VertexType, kMaxFaces,
    kMaxVertices>::setPlaneOfExistence(const Plane* a_plane) {
  plane_of_existence_m = a_plane;
}
template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
const Plane&
SegmentedHalfEdgePolygonCommon<FaceType, VertexType, kMaxFaces,
                               kMaxVertices>::getPlaneOfExistence(void) const {
  assert(plane_of_existence_m != nullptr);
  return *plane_of_existence_m;
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
Volume SegmentedHalfEdgePolygonCommon<FaceType, VertexType, kMaxFaces,
                                      kMaxVertices>::calculateVolume(void)
    const {
  return segmented_half_edge_polygon_detail::calculateMoments(
      this, Volume2D_Functor());
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
Volume SegmentedHalfEdgePolygonCommon<FaceType, VertexType, kMaxFaces,
                                      kMaxVertices>::calculateAbsoluteVolume(void)
    const {
  return std::fabs(this->calculateVolume());
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
Pt SegmentedHalfEdgePolygonCommon<FaceType, VertexType, kMaxFaces,
                                  kMaxVertices>::calculateCentroid(void) const {
  return segmented_half_edge_polygon_detail::calculateMoments(
      this, Centroid2D_Functor());
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
VolumeMoments
SegmentedHalfEdgePolygonCommon<FaceType, VertexType, kMaxFaces,
                               kMaxVertices>::calculateMoments(void) const {
  return segmented_half_edge_polygon_detail::calculateMoments(
      this, VolumeMoments2D_Functor());
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
VolumeMomentsAndNormal SegmentedHalfEdgePolygonCommon<
    FaceType, VertexType, kMaxFaces,
    kMaxVertices>::calculateVolumeMomentsAndNormal(void) const {
  return segmented_half_edge_polygon_detail::calculateMoments(
      this, VolumeMomentsAndNormal2D_Functor());
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
bool SegmentedHalfEdgePolygonCommon<
    FaceType, VertexType, kMaxFaces,
    kMaxVertices>::checkValidHalfEdgeStructure(void) {
  std::unordered_map<FaceType*, bool> face_found_through_opposite;
  for (auto& face : (*this)) {
    face_found_through_opposite[face] = false;
  }

  for (auto& face : (*this)) {
    bool valid_face = face->checkValidFace();
    if (!valid_face) {
      return false;
    }
  }

  if (!this->checkVerticesCompleteCycle()) {
    std::cout
        << "Vertex did not point to half edge that could cycle to itself. \n"
        << std::endl;
    return false;
  }

  if (!this->checkVerticesOnlyHaveHalfEdgesWithCorrectFaces()) {
    std::cout << "Contained a vertex that could cycle to a half-edge on a face"
              << " this polyhedron does not contain. \n"
              << std::endl;
    return false;
  }
  return true;
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
bool SegmentedHalfEdgePolygonCommon<
    FaceType, VertexType, kMaxFaces,
    kMaxVertices>::checkVerticesCompleteCycle(void) const {
  for (UnsignedIndex_t v = 0; v < this->getNumberOfVertices(); ++v) {
    if (!this->getVertex(v)->checkValidHalfEdgeCycle()) {
      return false;
    }
  }
  return true;
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
bool SegmentedHalfEdgePolygonCommon<
    FaceType, VertexType, kMaxFaces,
    kMaxVertices>::checkVerticesOnlyHaveHalfEdgesWithCorrectFaces(void) const {
  using HalfEdgeType = typename FaceType::half_edge_type;
  for (UnsignedIndex_t v = 0; v < this->getNumberOfVertices(); ++v) {
    const HalfEdgeType* current_half_edge = this->getVertex(v)->getHalfEdge();
    do {
      const auto current_half_edge_face = current_half_edge->getFace();
      if (current_half_edge_face == &getOpenBoundaryFace<FaceType>()) {
        // Check that the outside is still a valid loop.
        if (!FaceType::checkValidFaceLoop(&getOpenBoundaryFace<FaceType>(),
                                          current_half_edge)) {
          std::cout << "\n Polytope contains an OPEN_BOUNDARY_FACE that is "
                       "not a valid half-edge face. \n"
                    << std::endl;
          return false;
        }
      } else {
        // Make sure polytope contains this face
        auto iterator_of_face_in_polyhedron =
            std::find((*this).begin(), (*this).end(), current_half_edge_face);
        if (iterator_of_face_in_polyhedron == (*this).end() &&
            current_half_edge_face != &getOpenBoundaryFace<FaceType>()) {
          return false;
        }
      }
      current_half_edge =
          current_half_edge->getOppositeHalfEdge()->getPreviousHalfEdge();

    } while (current_half_edge != this->getVertex(v)->getHalfEdge());
  }
  return true;
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYGON_TPP_
