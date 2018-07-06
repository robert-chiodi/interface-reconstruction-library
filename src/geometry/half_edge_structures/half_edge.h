// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_H_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_H_

#include <float.h>
#include <utility>

#include "src/geometry/general/plane.h"
#include "src/geometry/general/pt.h"

namespace IRL {

// Forward declarations.
template <class PtType>
class Vertex;

template <class VertexType>
class HalfEdge;

template <class HalfEdgeType>
class Face;

template <class FaceType>
inline FaceType& getOpenBoundaryFace(void);

template <class VertexType>
inline void doubleLinkHalfEdges(HalfEdge<VertexType>* a_starting_half_edge,
                                HalfEdge<VertexType>* a_ending_half_edge);

template <class VertexType>
class HalfEdge {
 public:
  using vertex_type = VertexType;

  HalfEdge(void);

  HalfEdge(VertexType* a_vertex, HalfEdge* a_previous, HalfEdge* a_next,
           Face<HalfEdge>* a_face);

  HalfEdge(const HalfEdge& a_other) = default;
  HalfEdge& operator=(const HalfEdge& a_other) = default;

  void setPreviousHalfEdge(HalfEdge* a_previous);
  HalfEdge* getPreviousHalfEdge(void);
  const HalfEdge* getPreviousHalfEdge(void) const;

  void setNextHalfEdge(HalfEdge* a_next);
  HalfEdge* getNextHalfEdge(void);
  const HalfEdge* getNextHalfEdge(void) const;

  void setOppositeHalfEdge(HalfEdge* a_opposite);
  HalfEdge* getOppositeHalfEdge(void);
  const HalfEdge* getOppositeHalfEdge(void) const;

  void setFace(Face<HalfEdge>* a_face);
  Face<HalfEdge>* getFace(void);
  const Face<HalfEdge>* getFace(void) const;

  void setVertex(VertexType* a_vertex);
  VertexType* getVertex(void);
  const VertexType* getVertex(void) const;

  VertexType* getPreviousVertex(void);
  const VertexType* getPreviousVertex(void) const;

  ~HalfEdge(void) = default;

 private:
  HalfEdge* previous_m;
  HalfEdge* next_m;
  HalfEdge* opposite_m;
  VertexType* end_point_m;
  Face<HalfEdge>* corresponding_face_m;
};

template <class PtType>
class Vertex {
 public:
  using pt_type = PtType;

  Vertex(void);

  explicit Vertex(const PtType& a_location);

  void setHalfEdge(HalfEdge<Vertex>* a_half_edge);
  HalfEdge<Vertex>* getHalfEdge(void);
  const HalfEdge<Vertex>* getHalfEdge(void) const;

  const pt_type& getLocation(void) const;
  void setLocation(const PtType& a_location);

  void calculateDistanceToPlane(const Plane& a_plane);
  double getDistance(void) const;
  double setDistance(void) const;

  void markToBeClipped(void);
  void markToBeNotClipped(void);
  void setClip(const bool a_is_clipped);
  bool isClipped(void) const;
  bool isNotClipped(void) const;

  void setAsUnnecessaryToSeek(void);
  void setToSeek(void);
  bool needsToSeek(void) const;
  bool doesNotNeedToSeek(void) const;

  inline bool checkValidHalfEdgeCycle(void) const;

 private:
  PtType vertex_location_m;
  HalfEdge<Vertex>* half_edge_m;  // HalfEdge that ends at this vertex
  double distance_m;
  bool is_clipped_m;
  bool needs_to_seek_m;
};

template <class HalfEdgeType>
class Face {
 public:
  using half_edge_type = HalfEdgeType;

  Face(void);

  explicit Face(HalfEdgeType* a_starting_half_edge);

  void setStartingHalfEdge(HalfEdgeType* a_starting_half_edge);
  HalfEdgeType* getStartingHalfEdge(void);
  const HalfEdgeType* getStartingHalfEdge(void) const;
  void reverseFaceOrientation(void);

  void markAsVisited(void);
  void markAsNotVisited(void);
  bool hasBeenVisited(void) const;
  bool hasNotBeenVisited(void) const;

  inline bool checkValidFace(void);

  static inline bool checkValidFaceLoop(
      const Face* a_face, const HalfEdgeType* a_starting_half_edge);

  static inline bool checkValidCircuitOnFace(
      const HalfEdgeType* a_starting_half_edge);

  static inline bool checkCorrectReflectionOfFaceHalfEdgeOpposites(
      const HalfEdgeType* a_starting_half_edge);

  static inline bool checkAllHalfEdgesPointBackToFace(
      const Face* a_face, const HalfEdgeType* a_starting_half_edge);

 private:
  HalfEdgeType* starting_half_edge_m;
  bool has_been_visited_m;
};

}  // namespace IRL

#include "src/geometry/half_edge_structures/half_edge.tpp"

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_H_
