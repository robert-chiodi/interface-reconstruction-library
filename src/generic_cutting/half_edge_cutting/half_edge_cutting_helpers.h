// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_HELPERS_H_
#define SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_HELPERS_H_

#include "src/geometry/half_edge_structures/half_edge.h"
#include "src/helpers/geometric_cutting_helpers.h"

namespace IRL {


// Creates a vertex by interpolating the distances to find the
// zero location, and then subdivides the edge, creating
// new half edges and adjusting connectivity to
// maintain a valid half edge structure.
//
//
//               intersection half edge
// old_vertex  .---------------------->. old_vertex
//
//
//
//                new         old
//               ------->   -------->
// old_vertex  .----------x----------->. old_vertex
//               <-------   <--------
//                 old         new
//  
template <class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType>
void subdivideEdge(HalfEdgeType* a_half_edge_with_intersection,
                   SegmentedHalfEdgePolytopeType* a_polytope,
                   HalfEdgePolytopeType* a_complete_polytope);

template <class HalfEdgeType, class SegmentedHalfEdgePolygonType,
          class HalfEdgePolytopeType>
inline HalfEdgeType* separateIntersectedHalfEdge(
    HalfEdgeType* a_half_edge_with_intersection,
    SegmentedHalfEdgePolygonType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope);

template <class HalfEdgeType, class HalfEdgePolytopeType>
inline void createOppositeHalfEdgeFromIntersection(
    HalfEdgeType* a_half_edge_with_intersection,
    HalfEdgeType* a_newly_created_half_edge,
    HalfEdgePolytopeType* a_complete_polytope);

}  // namespace IRL

#include "src/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.tpp"

#endif  // SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_HELPERS_H_
