// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_H_
#define IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_H_

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/paraboloid_intersection/surface_output.h"
#include "irl/geometry/general/geometry_type_traits.h"
#include "irl/paraboloid_reconstruction/aligned_paraboloid.h"
#include "irl/paraboloid_reconstruction/ellipse.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/paraboloid_reconstruction/parametrized_surface.h"

namespace IRL {

template <class SegmentedHalfEdgePolyhedronType>
void integrateClipFaceCorrections(
    const SegmentedHalfEdgePolyhedronType& a_polytope,
    const AlignedParaboloid& a_paraboloid,
    const UnsignedIndex_t a_starting_face);

bool checkIfNewEdgeOnSameNappe(const AlignedParaboloid& a_paraboloid,
                               const Plane& a_face_plane, const Pt& a_pt_0,
                               const Pt& a_pt_1);

template <class PtType, class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType>
void placeSingleIntercept(const PtType& a_intersection_location,
                          HalfEdgeType* a_half_edge_with_intersection,
                          SegmentedHalfEdgePolytopeType* a_polytope,
                          HalfEdgePolytopeType* a_complete_polytope);

template <class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType, class VertexType>
void placeDoubleIntercept(const StackVector<std::pair<VertexType, double>, 2>&
                              a_intersection_location,
                          HalfEdgeType* a_half_edge_with_intersection,
                          SegmentedHalfEdgePolytopeType* a_polytope,
                          HalfEdgePolytopeType* a_complete_polytope);

template <class VertexType>
bool vertexBelow(const VertexType& a_pt, const AlignedParaboloid& a_paraboloid);

template <class VertexType>
void checkAndFindIntercepts(
    const AlignedParaboloid& a_paraboloid, const VertexType& a_pt_0,
    const VertexType& a_pt_1,
    StackVector<std::pair<VertexType, double>, 2>* a_intercepts);

bool isPtBeforeIntersectionWithEdge(const std::array<double, 2>& a_test_pt,
                                    const Pt& a_vertex_0, const Pt& a_vertex_1);

template <class HalfEdgeType>
bool needsWedgeCorrection(const AlignedParaboloid& a_paraboloid,
                          HalfEdgeType* a_start, HalfEdgeType* a_end);

template <class ReturnType, class HalfEdgeType,
          class SurfaceOutputType = NoSurfaceOutput>
ReturnType computeNewEdgeSegmentContribution(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_ref_pt,
    const HalfEdgeType a_entry_half_edge, const HalfEdgeType a_exit_half_edge,
    const bool skip_first = false, SurfaceOutputType* a_surface = NULL);

template <class ReturnType, class HalfEdgeType,
          class SurfaceOutputType = NoSurfaceOutput>
ReturnType computeContribution(const AlignedParaboloid& a_aligned_paraboloid,
                               const Pt& a_ref_pt,
                               const HalfEdgeType a_entry_half_edge,
                               const HalfEdgeType a_exit_half_edge,
                               const bool skip_first = false,
                               SurfaceOutputType* a_surface = NULL);

template <class ReturnType, class HalfEdgeType,
          class SurfaceOutputType = NoSurfaceOutput>
ReturnType orientAndApplyWedgeCorrection(const AlignedParaboloid& a_paraboloid,
                                         HalfEdgeType* a_start,
                                         HalfEdgeType* a_end,
                                         SurfaceOutputType* a_surface = NULL);

template <class ReturnType, class HalfEdgeType>
ReturnType orientAndApplyWedgeCorrectionHyperbolic(
    const AlignedParaboloid& a_paraboloid, HalfEdgeType* a_start,
    HalfEdgeType* a_end);

template <class HalfEdgeType>
bool ellipseContainedInFace(const Ellipse& a_ellipse,
                            HalfEdgeType* const a_half_edge);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
findFaceOnlyIntersections(SegmentedHalfEdgePolyhedronType* a_polytope,
                          HalfEdgePolytopeType* a_complete_polytope,
                          const AlignedParaboloid& a_aligned_paraboloid);

template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, void>
findEdgeIntersectionsElliptic(SegmentedHalfEdgePolyhedronType* a_polytope,
                              HalfEdgePolytopeType* a_complete_polytope,
                              const AlignedParaboloid& a_aligned_paraboloid);

template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, void>
triangulateAllNewFaces(SegmentedHalfEdgePolyhedronType* a_polytope,
                       HalfEdgePolytopeType* a_complete_polytope,
                       const AlignedParaboloid& a_aligned_paraboloid,
                       const UnsignedIndex_t a_face_start,
                       const UnsignedIndex_t a_face_end);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class SurfaceOutputType = NoSurfaceOutput>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
formParaboloidIntersectionBases(SegmentedHalfEdgePolyhedronType* a_polytope,
                                HalfEdgePolytopeType* a_complete_polytope,
                                const AlignedParaboloid& a_aligned_paraboloid,
                                const UnsignedIndex_t a_nudge_iter = 0,
                                SurfaceOutputType* a_surface = nullptr);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
formParaboloidIntersectionBasesElliptic(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AlignedParaboloid& a_aligned_paraboloid);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
formParaboloidIntersectionBasesHyperbolic(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AlignedParaboloid& a_aligned_paraboloid);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class SurfaceOutputType = NoSurfaceOutput>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithParaboloid(SegmentedHalfEdgePolyhedronType* a_polytope,
                                  HalfEdgePolytopeType* a_complete_polytope,
                                  const AlignedParaboloid& a_paraboloid,
                                  SurfaceOutputType* a_surface = NULL);
}  // namespace IRL

#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection.tpp"

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_H_
