// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_TPP_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_TPP_

namespace IRL {

template <class FaceType>
inline FaceType& getOpenBoundaryFace(void) {
  static FaceType OPEN_BOUNDARY_FACE;
  return OPEN_BOUNDARY_FACE;
}

template <class VertexType>
inline void doubleLinkHalfEdges(HalfEdge<VertexType>* a_starting_half_edge,
                                HalfEdge<VertexType>* a_ending_half_edge) {
  a_starting_half_edge->setNextHalfEdge(a_ending_half_edge);
  a_ending_half_edge->setPreviousHalfEdge(a_starting_half_edge);
}

template <class VertexType>
inline HalfEdge<VertexType>::HalfEdge(void)
    : previous_m(nullptr),
      next_m(nullptr),
      opposite_m(nullptr),
      end_point_m(nullptr),
      corresponding_face_m(nullptr) {}

template <class VertexType>
inline HalfEdge<VertexType>::HalfEdge(VertexType* a_vertex, HalfEdge* a_previous,
                               HalfEdge* a_next, Face<HalfEdge>* a_face)
    : previous_m(a_previous),
      next_m(a_next),
      opposite_m(nullptr),
      end_point_m(a_vertex),
      corresponding_face_m(a_face) {}

template <class VertexType>
inline void HalfEdge<VertexType>::setPreviousHalfEdge(HalfEdge* a_previous) {
  previous_m = a_previous;
}
template <class VertexType>
inline HalfEdge<VertexType>* HalfEdge<VertexType>::getPreviousHalfEdge(void) {
  return previous_m;
}
template <class VertexType>
inline const HalfEdge<VertexType>* HalfEdge<VertexType>::getPreviousHalfEdge(
    void) const {
  return previous_m;
}

template <class VertexType>
inline void HalfEdge<VertexType>::setNextHalfEdge(HalfEdge* a_next) {
  next_m = a_next;
}
template <class VertexType>
inline HalfEdge<VertexType>* HalfEdge<VertexType>::getNextHalfEdge(void) {
  return next_m;
}
template <class VertexType>
inline const HalfEdge<VertexType>* HalfEdge<VertexType>::getNextHalfEdge(void) const {
  return next_m;
}

template <class VertexType>
inline void HalfEdge<VertexType>::setOppositeHalfEdge(HalfEdge* a_opposite) {
  opposite_m = a_opposite;
}
template <class VertexType>
inline HalfEdge<VertexType>* HalfEdge<VertexType>::getOppositeHalfEdge(void) {
  return opposite_m;
}
template <class VertexType>
inline const HalfEdge<VertexType>* HalfEdge<VertexType>::getOppositeHalfEdge(
    void) const {
  return opposite_m;
}

template <class VertexType>
inline void HalfEdge<VertexType>::setFace(Face<HalfEdge<VertexType>>* a_face) {
  corresponding_face_m = a_face;
}
template <class VertexType>
inline Face<HalfEdge<VertexType>>* HalfEdge<VertexType>::getFace(void) {
  return corresponding_face_m;
}
template <class VertexType>
inline const Face<HalfEdge<VertexType>>* HalfEdge<VertexType>::getFace(void) const {
  return corresponding_face_m;
}

template <class VertexType>
inline void HalfEdge<VertexType>::setVertex(VertexType* a_vertex) {
  end_point_m = a_vertex;
}
template <class VertexType>
inline VertexType* HalfEdge<VertexType>::getVertex(void) {
  return end_point_m;
}
template <class VertexType>
inline const VertexType* HalfEdge<VertexType>::getVertex(void) const {
  return end_point_m;
}

template <class VertexType>
inline VertexType* HalfEdge<VertexType>::getPreviousVertex(void) {
  assert(previous_m != nullptr);
  return previous_m->getVertex();
}
template <class VertexType>
inline const VertexType* HalfEdge<VertexType>::getPreviousVertex(void) const {
  assert(previous_m != nullptr);
  return previous_m->getVertex();
}

template <class PtType>
inline Vertex<PtType>::Vertex(void)
    : vertex_location_m(Pt(0.0, 0.0, 0.0)),
      half_edge_m(nullptr),
      distance_m(DBL_MAX),
      is_clipped_m{false},
      needs_to_seek_m(false) {}

template <class PtType>
inline Vertex<PtType>::Vertex(const PtType& a_location)
    : vertex_location_m(a_location),
      half_edge_m(nullptr),
      distance_m(DBL_MAX),
      is_clipped_m{false},
      needs_to_seek_m(false) {}

template <class PtType>
inline void Vertex<PtType>::setHalfEdge(HalfEdge<Vertex>* a_half_edge) {
  half_edge_m = a_half_edge;
}
template <class PtType>
inline HalfEdge<Vertex<PtType>>* Vertex<PtType>::getHalfEdge(void) {
  return half_edge_m;
}
template <class PtType>
inline const HalfEdge<Vertex<PtType>>* Vertex<PtType>::getHalfEdge(void) const {
  return half_edge_m;
}

template <class PtType>
inline const PtType& Vertex<PtType>::getLocation(void) const {
  return vertex_location_m;
}

template <class PtType>
inline void Vertex<PtType>::setLocation(const PtType& a_location) {
  vertex_location_m = a_location;
}

template <class PtType>
inline void Vertex<PtType>::calculateDistanceToPlane(const Plane& a_plane) {
  distance_m = a_plane.signedDistanceToPoint(this->getLocation());
}
template <class PtType>
inline double Vertex<PtType>::getDistance(void) const {
  return distance_m;
}

template <class PtType>
inline void Vertex<PtType>::markToBeClipped(void) {
  is_clipped_m = true;
}
template <class PtType>
inline void Vertex<PtType>::markToBeNotClipped(void) {
  is_clipped_m = false;
}
template <class PtType>
inline void Vertex<PtType>::setClip(const bool a_is_clipped) {
  is_clipped_m = a_is_clipped;
}
	
template <class PtType>
inline bool Vertex<PtType>::isClipped(void) const {
  return is_clipped_m;
}
template <class PtType>
inline bool Vertex<PtType>::isNotClipped(void) const {
  return !is_clipped_m;
}
template <class PtType>
inline void Vertex<PtType>::setAsUnnecessaryToSeek(void) {
  needs_to_seek_m = false;
}
template <class PtType>
inline void Vertex<PtType>::setToSeek(void) {
  needs_to_seek_m = true;
}
template <class PtType>
inline bool Vertex<PtType>::needsToSeek(void) const {
  return needs_to_seek_m;
}
template <class PtType>
inline bool Vertex<PtType>::doesNotNeedToSeek(void) const {
  return !needs_to_seek_m;
}

template <class PtType>
bool Vertex<PtType>::checkValidHalfEdgeCycle(void) const {
  assert(half_edge_m != nullptr);
  const HalfEdge<Vertex>* current_half_edge = this->getHalfEdge();
  // To avoid hanging if enters another complete circuit
  // and can never return to start.
  constexpr UnsignedIndex_t max_loop_size = 1000;
  UnsignedIndex_t number_of_vertices = 0;
  do {
    if (number_of_vertices > max_loop_size) {
      std::cout
          << "Half edge for vertex did not cycle around vertex back to itself. "
          << std::endl;
      return false;
    }
    ++number_of_vertices;
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != this->getHalfEdge());
  return true;
}

template <class HalfEdgeType>
inline Face<HalfEdgeType>::Face(void)
    : starting_half_edge_m(nullptr), has_been_visited_m{false} {}

template <class HalfEdgeType>
inline Face<HalfEdgeType>::Face(HalfEdgeType* a_starting_half_edge)
    : starting_half_edge_m(a_starting_half_edge), has_been_visited_m{false} {}

template <class HalfEdgeType>
inline void Face<HalfEdgeType>::setStartingHalfEdge(
    HalfEdgeType* a_starting_half_edge) {
  starting_half_edge_m = a_starting_half_edge;
}
template <class HalfEdgeType>
inline HalfEdgeType* Face<HalfEdgeType>::getStartingHalfEdge(void) {
  assert(starting_half_edge_m != nullptr);
  return starting_half_edge_m;
}
template <class HalfEdgeType>
inline const HalfEdgeType* Face<HalfEdgeType>::getStartingHalfEdge(void) const {
  assert(starting_half_edge_m != nullptr);
  return starting_half_edge_m;
}

template <class HalfEdgeType>
void Face<HalfEdgeType>::reverseFaceOrientation(void) {
  HalfEdgeType* current_half_edge = this->getStartingHalfEdge();
  auto starting_vertex = current_half_edge->getVertex();
  do {
    current_half_edge->setVertex(current_half_edge->getPreviousVertex());
    HalfEdgeType* temporary_half_edge =
        current_half_edge->getPreviousHalfEdge();
    current_half_edge->setPreviousHalfEdge(
        current_half_edge->getNextHalfEdge());
    current_half_edge->setNextHalfEdge(temporary_half_edge);
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != this->getStartingHalfEdge());
  this->getStartingHalfEdge()->getPreviousHalfEdge()->setVertex(
      starting_vertex);
}

template <class HalfEdgeType>
inline void Face<HalfEdgeType>::markAsVisited(void) {
  has_been_visited_m = true;
}
template <class HalfEdgeType>
inline void Face<HalfEdgeType>::markAsNotVisited(void) {
  has_been_visited_m = false;
}
template <class HalfEdgeType>
inline bool Face<HalfEdgeType>::hasBeenVisited(void) const {
  return has_been_visited_m;
}
template <class HalfEdgeType>
inline bool Face<HalfEdgeType>::hasNotBeenVisited(void) const {
  return !this->hasBeenVisited();
}

template <class HalfEdgeType>
inline bool Face<HalfEdgeType>::checkValidFace(void) {
  return Face::checkValidFaceLoop(this, this->getStartingHalfEdge());
}

template <class HalfEdgeType>
bool Face<HalfEdgeType>::checkValidFaceLoop(
    const Face* a_face, const HalfEdgeType* a_starting_half_edge) {
  return Face::checkValidCircuitOnFace(a_starting_half_edge) &&
         Face::checkCorrectReflectionOfFaceHalfEdgeOpposites(
             a_starting_half_edge) &&
         Face::checkAllHalfEdgesPointBackToFace(a_face, a_starting_half_edge);
}

template <class HalfEdgeType>
bool Face<HalfEdgeType>::checkValidCircuitOnFace(
    const HalfEdgeType* a_starting_half_edge) {
  assert(a_starting_half_edge != nullptr);
  const HalfEdgeType* current_half_edge = a_starting_half_edge;
  // To avoid hanging if enters another complete circuit
  // and can never return to start.
  constexpr UnsignedIndex_t max_loop_size = 1000;
  UnsignedIndex_t number_of_vertices = 0;
  do {
    const HalfEdgeType* next_half_edge = current_half_edge->getNextHalfEdge();
    if (next_half_edge == nullptr) {
      std::cout << "Nullptr present in HalfEdge chain on face. " << std::endl;
      return false;
    }
    if (next_half_edge->getPreviousHalfEdge() != current_half_edge) {
      std::cout << "Half edge goes from "
                << current_half_edge->getPreviousVertex()->getLocation()
                << " to " << current_half_edge->getVertex()->getLocation()
                << std::endl;
      std::cout << "Incorrect double linking between successive half edges."
                << std::endl;
      std::cout << "Current half edge is at location "
                << current_half_edge->getVertex()->getLocation()
                << " with half edge address " << current_half_edge << std::endl;
      std::cout << "Next half edge has previous location as "
                << next_half_edge->getPreviousVertex()->getLocation()
                << " with previous half edge address of "
                << next_half_edge->getPreviousHalfEdge() << std::endl;
      std::cout << "Two vertex addresses are " << current_half_edge->getVertex()
                << " and " << next_half_edge->getPreviousVertex()
                << " , repsectively. " << std::endl;
      return false;
    }
    if (number_of_vertices > max_loop_size) {
      std::cout << "Apparent incorrect circuit of HalfEdges on face. "
                << std::endl;
      return false;
    }
    ++number_of_vertices;
    current_half_edge = next_half_edge;
  } while (current_half_edge != a_starting_half_edge);
  return true;
}

template <class HalfEdgeType>
bool Face<HalfEdgeType>::checkCorrectReflectionOfFaceHalfEdgeOpposites(
    const HalfEdgeType* a_starting_half_edge) {
  assert(a_starting_half_edge != nullptr);
  const HalfEdgeType* current_half_edge = a_starting_half_edge;
  // To avoid hanging if enters another complete circuit
  // and can never return to start.
  constexpr UnsignedIndex_t max_loop_size = 10000;
  UnsignedIndex_t number_of_vertices = 0;
  do {
    assert(current_half_edge != nullptr);
    assert(current_half_edge->getOppositeHalfEdge() != nullptr);
    const HalfEdgeType* opposite_half_edge =
        current_half_edge->getOppositeHalfEdge();
    if (opposite_half_edge->getOppositeHalfEdge() != current_half_edge) {
      std::cout
          << "Half edge in face is not the opposite of its opposite half edge. "
          << std::endl;
      return false;
    }

    if (number_of_vertices > max_loop_size) {
      std::cout << "Apparent incorrect circuit of HalfEdges on face. "
                << std::endl;
      return false;
    }
    ++number_of_vertices;
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != a_starting_half_edge);
  return true;
}

template <class HalfEdgeType>
bool Face<HalfEdgeType>::checkAllHalfEdgesPointBackToFace(
    const Face* a_face, const HalfEdgeType* a_starting_half_edge) {
  assert(a_starting_half_edge != nullptr);
  const HalfEdgeType* current_half_edge = a_starting_half_edge;
  // To avoid hanging if enters another complete circuit
  // and can never return to start.
  constexpr UnsignedIndex_t max_loop_size = 10000;
  UnsignedIndex_t number_of_vertices = 0;
  do {
    if (current_half_edge->getFace() != a_face) {
      std::cout << current_half_edge->getFace() << std::endl;
      std::cout << current_half_edge->getVertex()->getLocation() << std::endl;
      std::cout << "Half edge does not point to the face it is actually on. "
                << std::endl;
      return false;
    }

    if (number_of_vertices > max_loop_size) {
      std::cout << "Apparent incorrect circuit of HalfEdges on face. "
                << std::endl;
      return false;
    }
    ++number_of_vertices;
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != a_starting_half_edge);
  return true;
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_TPP_
