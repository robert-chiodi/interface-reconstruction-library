// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// NOTE: THIS SHOULD BE WRITTEN WITH JUST THE OTHER HALF EDGE, WHERE
// A NEW FACE CLASS IS MADE AND EVERYTHING TEMPLATED ON IT. THIS
// HOWEVER LEADS TO A RECURSIVE INSTANTIATION OF TEMPLATES, SINCE
// THE HALF EDGE CLASS AND FACE CLASS RELY ON ONE ANOTHER.
// FOR NOW, JUST HARD CODE DIFFERENT CLASSES.

#ifndef IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_PARABOLOID_H_
#define IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_PARABOLOID_H_

#include <float.h>
#include <utility>

#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/half_edge_structures/half_edge.h"

namespace IRL {

// Forward declarations.
template <class PtType>
class VertexParaboloid;

template <class VertexType>
class HalfEdgeParaboloid;

template <class HalfEdgeType>
class FaceParaboloid;

template <class VertexType>
inline void doubleLinkHalfEdges(
    HalfEdgeParaboloid<VertexType>* a_starting_half_edge,
    HalfEdgeParaboloid<VertexType>* a_ending_half_edge);

template <class VertexType>
inline void setMutualOpposites(
    HalfEdgeParaboloid<VertexType>* a_starting_half_edge,
    HalfEdgeParaboloid<VertexType>* a_ending_half_edge);

template <class VertexType>
class HalfEdgeParaboloid {
 public:
  using vertex_type = VertexType;
  using value_type = typename vertex_type::value_type;

  HalfEdgeParaboloid(void);

  HalfEdgeParaboloid(VertexType* a_vertex, HalfEdgeParaboloid* a_previous,
                     HalfEdgeParaboloid* a_next,
                     FaceParaboloid<HalfEdgeParaboloid>* a_face);

  HalfEdgeParaboloid(const HalfEdgeParaboloid& a_other) = default;
  HalfEdgeParaboloid& operator=(const HalfEdgeParaboloid& a_other) = default;

  void setPreviousHalfEdge(HalfEdgeParaboloid* a_previous);
  HalfEdgeParaboloid* getPreviousHalfEdge(void);
  const HalfEdgeParaboloid* getPreviousHalfEdge(void) const;

  void setNextHalfEdge(HalfEdgeParaboloid* a_next);
  HalfEdgeParaboloid* getNextHalfEdge(void);
  const HalfEdgeParaboloid* getNextHalfEdge(void) const;

  void setOppositeHalfEdge(HalfEdgeParaboloid* a_opposite);
  HalfEdgeParaboloid* getOppositeHalfEdge(void);
  const HalfEdgeParaboloid* getOppositeHalfEdge(void) const;

  void setFace(FaceParaboloid<HalfEdgeParaboloid>* a_face);
  FaceParaboloid<HalfEdgeParaboloid>* getFace(void);
  const FaceParaboloid<HalfEdgeParaboloid>* getFace(void) const;

  void setVertex(VertexType* a_vertex);
  VertexType* getVertex(void);
  const VertexType* getVertex(void) const;

  VertexType* getPreviousVertex(void);
  const VertexType* getPreviousVertex(void) const;

  ~HalfEdgeParaboloid(void) = default;

 private:
  HalfEdgeParaboloid* previous_m;
  HalfEdgeParaboloid* next_m;
  HalfEdgeParaboloid* opposite_m;
  VertexType* end_point_m;
  FaceParaboloid<HalfEdgeParaboloid>* corresponding_face_m;
};

template <class PtType>
class VertexParaboloid {
 public:
  using pt_type = PtType;
  using value_type = typename pt_type::value_type;

  VertexParaboloid(void);

  explicit VertexParaboloid(const PtType& a_location);

  void setHalfEdge(HalfEdgeParaboloid<VertexParaboloid>* a_half_edge);
  HalfEdgeParaboloid<VertexParaboloid>* getHalfEdge(void);
  const HalfEdgeParaboloid<VertexParaboloid>* getHalfEdge(void) const;

  pt_type& getLocation(void);
  const pt_type& getLocation(void) const;
  void setLocation(const PtType& a_location);

  void calculateDistanceToPlane(const PlaneBase<value_type>& a_plane);
  value_type getDistance(void) const;
  value_type setDistance(void) const;

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
  HalfEdgeParaboloid<VertexParaboloid>*
      half_edge_m;  // HalfEdgeParaboloid that ends at this vertex
  value_type distance_m;
  bool is_clipped_m = false;
  bool needs_to_seek_m = false;
};

template <class HalfEdgeType>
class FaceParaboloid : public Face<HalfEdgeType> {
 public:
  using half_edge_type = HalfEdgeType;
  using value_type = typename half_edge_type::value_type;

  FaceParaboloid(void);

  explicit FaceParaboloid(HalfEdgeType* a_starting_half_edge);

  void setPlane(const PlaneBase<value_type>& a_plane);
  const PlaneBase<value_type>& getPlane(void) const;

  void clearIntersections(void);
  void addIntersection(void);
  void addDoubleIntersection(void);
  UnsignedIndex_t getNumberOfIntersections(void) const;

  /// \brief This is a subset of the intersections whose edge is directly
  /// parallel to the tangent of the paraboloid at that point.
  void addEdgeParallelIntersection(void);
  void addEdgeParallelIntersections(const UnsignedIndex_t a_intersections);
  UnsignedIndex_t getNumberOfEdgeParallelIntersections(void) const;
  void setAsTriangle(void);
  void setAsNotTriangle(void);
  bool isTriangle(void) const;

 private:
  PlaneBase<value_type> face_plane_m;
  UnsignedIndex_t intersections_m;
  UnsignedIndex_t edge_parallel_intersections_m;
  bool is_triangle_m = false;
};

}  // namespace IRL

#include "irl/geometry/half_edge_structures/half_edge_paraboloid.tpp"

#endif  // IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_PARABOLOID_H_
