// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_SIMPLEX_CUTTING_SIMPLEX_CUTTING_H_
#define SRC_GENERIC_CUTTING_SIMPLEX_CUTTING_SIMPLEX_CUTTING_H_

#include <algorithm>
#include <utility>

#include "src/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "src/geometry/general/geometry_type_traits.h"
#include "src/geometry/half_edge_structures/half_edge.h"
#include "src/helpers/SFINAE_boiler_plate.h"

namespace IRL {

// The idea here is to take a complete_polytope, and partition it amongst
// reconstructions to obtain a_polytope, which is a subset of
// a_complete_polytope. The HalfEdgePolytope a_complete_polytope provides
// the underlying data, where new HalfEdge, Pt, Vertex, and Face objects
// should be obtained from.
template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value>
splitDecomposedPolytope(SegmentedHalfEdgePolyhedronType* a_polytope,
                        SegmentedHalfEdgePolyhedronType* a_clipped_polytope,
                        HalfEdgePolytopeType* a_complete_polytope,
                        const Plane& a_plane);

template <class SegmentedHalfEdgePolygonType, class HalfEdgePolytopeType>
enable_if_t<is_polygon<SegmentedHalfEdgePolygonType>::value>
splitDecomposedPolytope(SegmentedHalfEdgePolygonType* a_polytope,
                        SegmentedHalfEdgePolygonType* a_clipped_polytope,
                        HalfEdgePolytopeType* a_complete_polytope,
                        const Plane& a_plane);

template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value>
truncateDecomposedPolytope(SegmentedHalfEdgePolyhedronType* a_polytope,
                           HalfEdgePolytopeType* a_complete_polytope,
                           const Plane& a_plane);

template <class SegmentedHalfEdgePolygonType, class HalfEdgePolytopeType>
enable_if_t<is_polygon<SegmentedHalfEdgePolygonType>::value>
truncateDecomposedPolytope(SegmentedHalfEdgePolygonType* a_polytope,
                           HalfEdgePolytopeType* a_complete_polytope,
                           const Plane& a_plane);

}  // namespace IRL

#include "src/generic_cutting/simplex_cutting/simplex_cutting.tpp"

#endif  // SRC_GENERIC_CUTTING_SIMPLEX_CUTTING_SIMPLEX_CUTTING_H_
