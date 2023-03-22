// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_HELPERS_TPP_
#define IRL_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_HELPERS_TPP_

namespace IRL {

template <class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType>
void subdivideEdge(HalfEdgeType* a_half_edge_with_intersection,
                   SegmentedHalfEdgePolytopeType* a_polytope,
                   HalfEdgePolytopeType* a_complete_polytope) {
  const auto previous_half_edge =
      a_half_edge_with_intersection->getPreviousHalfEdge();
  auto first_intersection_vertex = a_complete_polytope->getNewVertex(
      typename SegmentedHalfEdgePolytopeType::vertex_type(
          a_half_edge_with_intersection->getVertex()
              ->getLocation()
              .fromEdgeIntersection(
                  previous_half_edge->getVertex()->getLocation(),
                  previous_half_edge->getVertex()->getDistance(),
                  a_half_edge_with_intersection->getVertex()->getLocation(),
                  a_half_edge_with_intersection->getVertex()->getDistance())));
  a_polytope->addVertex(first_intersection_vertex);

  HalfEdgeType* new_half_edge = separateIntersectedHalfEdge(
      first_intersection_vertex, a_half_edge_with_intersection, a_polytope,
      a_complete_polytope);
  createOppositeHalfEdgeFromIntersection(a_half_edge_with_intersection,
                                         new_half_edge, a_complete_polytope);
}

template <class VertexType, class HalfEdgeType,
          class SegmentedHalfEdgePolytopeType, class HalfEdgePolytopeType>
inline HalfEdgeType* separateIntersectedHalfEdge(
    VertexType* a_intersection_vertex,
    HalfEdgeType* a_half_edge_with_intersection,
    SegmentedHalfEdgePolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope) {
  HalfEdgeType* previous_half_edge =
      a_half_edge_with_intersection->getPreviousHalfEdge();

  HalfEdgeType* first_intersection_half_edge =
      a_complete_polytope->getNewHalfEdge(
          HalfEdgeType(a_intersection_vertex, previous_half_edge,
                       a_half_edge_with_intersection,
                       a_half_edge_with_intersection->getFace()));

  a_intersection_vertex->setHalfEdge(first_intersection_half_edge);

  previous_half_edge->setNextHalfEdge(first_intersection_half_edge);
  a_half_edge_with_intersection->setPreviousHalfEdge(
      first_intersection_half_edge);

  // Note: This does not set the Opposite HalfEdge yet!

  return first_intersection_half_edge;
}

template <class HalfEdgeType, class HalfEdgePolytopeType>
inline void createOppositeHalfEdgeFromIntersection(
    HalfEdgeType* a_half_edge_with_intersection,
    HalfEdgeType* a_newly_created_half_edge,
    HalfEdgePolytopeType* a_complete_polytope) {
  HalfEdgeType* opposite_half_edge_of_intersection =
      a_half_edge_with_intersection->getOppositeHalfEdge();

  HalfEdgeType* opposite_new_half_edge_from_intersection =
      a_complete_polytope->getNewHalfEdge(HalfEdgeType(
          a_newly_created_half_edge->getVertex(),
          opposite_half_edge_of_intersection->getPreviousHalfEdge(),
          opposite_half_edge_of_intersection,
          opposite_half_edge_of_intersection->getFace()));

  opposite_half_edge_of_intersection->getPreviousHalfEdge()->setNextHalfEdge(
      opposite_new_half_edge_from_intersection);
  opposite_half_edge_of_intersection->setPreviousHalfEdge(
      opposite_new_half_edge_from_intersection);

  // Set opposites now
  a_half_edge_with_intersection->setOppositeHalfEdge(
      opposite_new_half_edge_from_intersection);
  opposite_new_half_edge_from_intersection->setOppositeHalfEdge(
      a_half_edge_with_intersection);

  a_newly_created_half_edge->setOppositeHalfEdge(
      opposite_half_edge_of_intersection);
  opposite_half_edge_of_intersection->setOppositeHalfEdge(
      a_newly_created_half_edge);
}

}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_HELPERS_TPP_
