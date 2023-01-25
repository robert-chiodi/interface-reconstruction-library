// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_PARABOLOID_TPP_
#define IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_PARABOLOID_TPP_

namespace IRL {

template <class VertexType>
inline void doubleLinkHalfEdges(
    HalfEdgeParaboloid<VertexType>* a_starting_half_edge,
    HalfEdgeParaboloid<VertexType>* a_ending_half_edge) {
  a_starting_half_edge->setNextHalfEdge(a_ending_half_edge);
  a_ending_half_edge->setPreviousHalfEdge(a_starting_half_edge);
}

template <class VertexType>
inline void setMutualOpposites(HalfEdgeParaboloid<VertexType>* a_half_edge_0,
                               HalfEdgeParaboloid<VertexType>* a_half_edge_1) {
  a_half_edge_0->setOppositeHalfEdge(a_half_edge_1);
  a_half_edge_1->setOppositeHalfEdge(a_half_edge_0);
}

template <class VertexType>
inline HalfEdgeParaboloid<VertexType>::HalfEdgeParaboloid(void)
    : previous_m(nullptr),
      next_m(nullptr),
      opposite_m(nullptr),
      end_point_m(nullptr),
      corresponding_face_m(nullptr) {}

template <class VertexType>
inline HalfEdgeParaboloid<VertexType>::HalfEdgeParaboloid(
    VertexType* a_vertex, HalfEdgeParaboloid* a_previous,
    HalfEdgeParaboloid* a_next, FaceParaboloid<HalfEdgeParaboloid>* a_face)
    : previous_m(a_previous),
      next_m(a_next),
      opposite_m(nullptr),
      end_point_m(a_vertex),
      corresponding_face_m(a_face) {}

template <class VertexType>
inline void HalfEdgeParaboloid<VertexType>::setPreviousHalfEdge(
    HalfEdgeParaboloid* a_previous) {
  previous_m = a_previous;
}
template <class VertexType>
inline HalfEdgeParaboloid<VertexType>*
HalfEdgeParaboloid<VertexType>::getPreviousHalfEdge(void) {
  return previous_m;
}
template <class VertexType>
inline const HalfEdgeParaboloid<VertexType>*
HalfEdgeParaboloid<VertexType>::getPreviousHalfEdge(void) const {
  return previous_m;
}

template <class VertexType>
inline void HalfEdgeParaboloid<VertexType>::setNextHalfEdge(
    HalfEdgeParaboloid* a_next) {
  next_m = a_next;
}
template <class VertexType>
inline HalfEdgeParaboloid<VertexType>*
HalfEdgeParaboloid<VertexType>::getNextHalfEdge(void) {
  return next_m;
}
template <class VertexType>
inline const HalfEdgeParaboloid<VertexType>*
HalfEdgeParaboloid<VertexType>::getNextHalfEdge(void) const {
  return next_m;
}

template <class VertexType>
inline void HalfEdgeParaboloid<VertexType>::setOppositeHalfEdge(
    HalfEdgeParaboloid* a_opposite) {
  opposite_m = a_opposite;
}
template <class VertexType>
inline HalfEdgeParaboloid<VertexType>*
HalfEdgeParaboloid<VertexType>::getOppositeHalfEdge(void) {
  return opposite_m;
}
template <class VertexType>
inline const HalfEdgeParaboloid<VertexType>*
HalfEdgeParaboloid<VertexType>::getOppositeHalfEdge(void) const {
  return opposite_m;
}

template <class VertexType>
inline void HalfEdgeParaboloid<VertexType>::setFace(
    FaceParaboloid<HalfEdgeParaboloid<VertexType>>* a_face) {
  corresponding_face_m = a_face;
}
template <class VertexType>
inline FaceParaboloid<HalfEdgeParaboloid<VertexType>>*
HalfEdgeParaboloid<VertexType>::getFace(void) {
  return corresponding_face_m;
}
template <class VertexType>
inline const FaceParaboloid<HalfEdgeParaboloid<VertexType>>*
HalfEdgeParaboloid<VertexType>::getFace(void) const {
  return corresponding_face_m;
}

template <class VertexType>
inline void HalfEdgeParaboloid<VertexType>::setVertex(VertexType* a_vertex) {
  end_point_m = a_vertex;
}
template <class VertexType>
inline VertexType* HalfEdgeParaboloid<VertexType>::getVertex(void) {
  return end_point_m;
}
template <class VertexType>
inline const VertexType* HalfEdgeParaboloid<VertexType>::getVertex(void) const {
  return end_point_m;
}

template <class VertexType>
inline VertexType* HalfEdgeParaboloid<VertexType>::getPreviousVertex(void) {
  assert(previous_m != nullptr);
  return previous_m->getVertex();
}
template <class VertexType>
inline const VertexType* HalfEdgeParaboloid<VertexType>::getPreviousVertex(
    void) const {
  assert(previous_m != nullptr);
  return previous_m->getVertex();
}

template <class PtType>
inline VertexParaboloid<PtType>::VertexParaboloid(void)
    : vertex_location_m(PtBase<typename PtType::value_type>(
          static_cast<typename PtType::value_type>(0),
          static_cast<typename PtType::value_type>(0),
          static_cast<typename PtType::value_type>(0))),
      half_edge_m(nullptr),
      distance_m(DBL_MAX),
      is_clipped_m{false},
      needs_to_seek_m(false),
      is_entry_m(false) {}

template <class PtType>
inline VertexParaboloid<PtType>::VertexParaboloid(const PtType& a_location)
    : vertex_location_m(a_location),
      half_edge_m(nullptr),
      distance_m(DBL_MAX),
      is_clipped_m{false},
      needs_to_seek_m(false),
      is_entry_m(false) {}

template <class PtType>
inline void VertexParaboloid<PtType>::setHalfEdge(
    HalfEdgeParaboloid<VertexParaboloid>* a_half_edge) {
  half_edge_m = a_half_edge;
}
template <class PtType>
inline HalfEdgeParaboloid<VertexParaboloid<PtType>>*
VertexParaboloid<PtType>::getHalfEdge(void) {
  return half_edge_m;
}
template <class PtType>
inline const HalfEdgeParaboloid<VertexParaboloid<PtType>>*
VertexParaboloid<PtType>::getHalfEdge(void) const {
  return half_edge_m;
}

template <class PtType>
inline PtType& VertexParaboloid<PtType>::getLocation(void) {
  return vertex_location_m;
}

template <class PtType>
inline const PtType& VertexParaboloid<PtType>::getLocation(void) const {
  return vertex_location_m;
}

template <class PtType>
inline void VertexParaboloid<PtType>::setLocation(const PtType& a_location) {
  vertex_location_m = a_location;
}

template <class PtType>
inline void VertexParaboloid<PtType>::calculateDistanceToPlane(
    const PlaneBase<typename PtType::value_type>& a_plane) {
  distance_m = a_plane.signedDistanceToPoint(this->getLocation());
}
template <class PtType>
inline typename PtType::value_type VertexParaboloid<PtType>::getDistance(
    void) const {
  return distance_m;
}

template <class PtType>
inline void VertexParaboloid<PtType>::markToBeClipped(void) {
  is_clipped_m = true;
}
template <class PtType>
inline void VertexParaboloid<PtType>::markToBeNotClipped(void) {
  is_clipped_m = false;
}
template <class PtType>
inline void VertexParaboloid<PtType>::setClip(const bool a_is_clipped) {
  is_clipped_m = a_is_clipped;
}

template <class PtType>
inline bool VertexParaboloid<PtType>::isClipped(void) const {
  return is_clipped_m;
}
template <class PtType>
inline bool VertexParaboloid<PtType>::isNotClipped(void) const {
  return !is_clipped_m;
}
template <class PtType>
inline void VertexParaboloid<PtType>::setAsUnnecessaryToSeek(void) {
  needs_to_seek_m = false;
}
template <class PtType>
inline void VertexParaboloid<PtType>::setToSeek(void) {
  needs_to_seek_m = true;
}
template <class PtType>
inline bool VertexParaboloid<PtType>::needsToSeek(void) const {
  return needs_to_seek_m;
}
template <class PtType>
inline bool VertexParaboloid<PtType>::doesNotNeedToSeek(void) const {
  return !needs_to_seek_m;
}
template <class PtType>
inline void VertexParaboloid<PtType>::markAsEntry(void) {
  is_entry_m = true;
}
template <class PtType>
inline void VertexParaboloid<PtType>::markAsExit(void) {
  is_entry_m = false;
}
template <class PtType>
inline bool VertexParaboloid<PtType>::isEntry(void) const {
  return is_entry_m;
}
template <class PtType>
inline bool VertexParaboloid<PtType>::isExit(void) const {
  return !is_entry_m;
}

template <class PtType>
bool VertexParaboloid<PtType>::checkValidHalfEdgeCycle(void) const {
  assert(half_edge_m != nullptr);
  const HalfEdgeParaboloid<VertexParaboloid>* current_half_edge =
      this->getHalfEdge();
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
inline FaceParaboloid<HalfEdgeType>::FaceParaboloid(void)
    : Face<HalfEdgeType>(),
      face_plane_m(),
      intersections_m(0),
      edge_parallel_intersections_m(0) {}

template <class HalfEdgeType>
inline FaceParaboloid<HalfEdgeType>::FaceParaboloid(
    HalfEdgeType* a_starting_half_edge)
    : Face<HalfEdgeType>(a_starting_half_edge),
      face_plane_m(),
      intersections_m(0),
      edge_parallel_intersections_m(0) {}

template <class HalfEdgeType>
inline void FaceParaboloid<HalfEdgeType>::setPlane(
    const PlaneBase<typename HalfEdgeType::value_type>& a_plane) {
  face_plane_m = a_plane;
}
template <class HalfEdgeType>
inline const PlaneBase<typename HalfEdgeType::value_type>&
FaceParaboloid<HalfEdgeType>::getPlane(void) const {
  return face_plane_m;
}

template <class HalfEdgeType>
void FaceParaboloid<HalfEdgeType>::clearIntersections(void) {
  intersections_m = 0;
  edge_parallel_intersections_m = 0;
}

template <class HalfEdgeType>
void FaceParaboloid<HalfEdgeType>::addIntersection(void) {
  ++intersections_m;
}

template <class HalfEdgeType>
void FaceParaboloid<HalfEdgeType>::addDoubleIntersection(void) {
  intersections_m += 2;
}

template <class HalfEdgeType>
UnsignedIndex_t FaceParaboloid<HalfEdgeType>::getNumberOfIntersections(
    void) const {
  return intersections_m;
}

template <class HalfEdgeType>
void FaceParaboloid<HalfEdgeType>::addEdgeParallelIntersection(void) {
  ++edge_parallel_intersections_m;
}

template <class HalfEdgeType>
void FaceParaboloid<HalfEdgeType>::addEdgeParallelIntersections(
    const UnsignedIndex_t a_intersections) {
  edge_parallel_intersections_m += a_intersections;
}

template <class HalfEdgeType>
UnsignedIndex_t
FaceParaboloid<HalfEdgeType>::getNumberOfEdgeParallelIntersections(void) const {
  return edge_parallel_intersections_m;
}

}  // namespace IRL

#endif  // IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_PARABOLOID_TPP_
