// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYHEDRON_TPP_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYHEDRON_TPP_

#include "src/geometry/general/moment_calculation_through_simplices.h"

namespace IRL {

namespace segmented_half_edge_polyhedron_detail {
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
  if (a_geometry->getNumberOfFaces() == 0) {
    return a_moment_accumulator.getMoments();
  }

  TetPtReferenceWrapper<typename GeometryType::pt_type> facet_and_pt;
  // Mark all faces as unvisited
  for (auto& face : (*a_geometry)) {
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

}  // namespace segmented_half_edge_polyhedron_detail

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
Volume SegmentedHalfEdgePolyhedronCommon<FaceType, VertexType, kMaxFaces,
                                         kMaxVertices>::calculateVolume(void) {
  return segmented_half_edge_polyhedron_detail::calculateMoments(
      this, Volume3D_Functor());
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
Volume
SegmentedHalfEdgePolyhedronCommon<FaceType, VertexType, kMaxFaces,
                                  kMaxVertices>::calculateAbsoluteVolume(void) {
  return std::fabs(this->calculateVolume());
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
Pt SegmentedHalfEdgePolyhedronCommon<FaceType, VertexType, kMaxFaces,
                                     kMaxVertices>::calculateCentroid(void) {
  return segmented_half_edge_polyhedron_detail::calculateMoments(
      this, Centroid3D_Functor());
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
VolumeMoments SegmentedHalfEdgePolyhedronCommon<
    FaceType, VertexType, kMaxFaces, kMaxVertices>::calculateMoments(void) {
  return segmented_half_edge_polyhedron_detail::calculateMoments(
      this, VolumeMoments3D_Functor());
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
bool SegmentedHalfEdgePolyhedronCommon<
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
    this->markFacesTouchedByOpposite(face, &face_found_through_opposite);
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

  for (const auto& element : face_found_through_opposite) {
    if (!element.second) {
      std::cout << "Face exists in polyhedron that was not accessed through "
                   "an opposite of a half-edge"
                << std::endl;
      return false;
    }
  }
  return true;
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
void SegmentedHalfEdgePolyhedronCommon<FaceType, VertexType, kMaxFaces,
                                       kMaxVertices>::
    markFacesTouchedByOpposite(
        FaceType* a_face_to_touch_from,
        std::unordered_map<FaceType*, bool>* a_map_to_fill) {
  auto current_half_edge = a_face_to_touch_from->getStartingHalfEdge();
  do {
    (*a_map_to_fill)[current_half_edge->getOppositeHalfEdge()->getFace()] =
        true;
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != a_face_to_touch_from->getStartingHalfEdge());
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
bool SegmentedHalfEdgePolyhedronCommon<
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
bool SegmentedHalfEdgePolyhedronCommon<
    FaceType, VertexType, kMaxFaces,
    kMaxVertices>::checkVerticesOnlyHaveHalfEdgesWithCorrectFaces(void) const {
  using HalfEdgeType = typename FaceType::half_edge_type;
  for (UnsignedIndex_t v = 0; v < this->getNumberOfVertices(); ++v) {
    const HalfEdgeType* current_half_edge = this->getVertex(v)->getHalfEdge();
    do {
      const auto current_half_edge_face = current_half_edge->getFace();
      auto iterator_of_face_in_polyhedron =
          std::find((*this).begin(), (*this).end(), current_half_edge_face);
      if (iterator_of_face_in_polyhedron == (*this).end()) {
        std::cout << "Face not belonging: " << current_half_edge_face
                  << std::endl;
        std::cout << "Faces in polyhedron: " << std::endl;
        for (const auto& face : (*this)) {
          std::cout << face << std::endl;
        }
        std::cout << "Vertex being cycle " << this->getVertex(v)->getLocation()
                  << std::endl;
        std::cout << "Half edges on this vertex come from: " << std::endl;
        const HalfEdgeType* error_half_edge_cycling = current_half_edge;
        do {
          std::cout
              << error_half_edge_cycling->getPreviousVertex()->getLocation()
              << std::endl;
          error_half_edge_cycling =
              error_half_edge_cycling->getOppositeHalfEdge()
                  ->getPreviousHalfEdge();
        } while (error_half_edge_cycling != current_half_edge);
        return false;
      }
      current_half_edge =
          current_half_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
    } while (current_half_edge != this->getVertex(v)->getHalfEdge());
  }
  return true;
}

template <class FunctorType, UnsignedIndex_t kArrayLength, class FaceType,
          class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
inline VolumeMomentsAndDoubles<kArrayLength>
SegmentedHalfEdgePolyhedronSpecificPt<
    PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>, FaceType,
    VertexType, kMaxFaces,
    kMaxVertices>::calculateVolumeMomentsAndDoubles(void) {
  return segmented_half_edge_polyhedron_detail::calculateMoments(
      this, VolumeMomentsAndDoubles3D_Functor<kArrayLength>());
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYHEDRON_TPP_
