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

#include <float.h>
#include <cassert>
#include <cmath>

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/generic_cutting/paraboloid_intersection/surface_output.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/helpers/mymath.h"
#include "irl/moments/volume_moments_with_gradient.h"
#include "irl/moments/volume_with_gradient.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/paraboloid_reconstruction/parametrized_surface.h"
#include "irl/paraboloid_reconstruction/rational_bezier_arc.h"

namespace IRL {

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class ParaboloidType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithParaboloid(SegmentedHalfEdgePolyhedronType* a_polytope,
                                  HalfEdgePolytopeType* a_complete_polytope,
                                  const ParaboloidType& a_paraboloid);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class AlignedParaboloidType,
          class SurfaceOutputType = NoSurfaceOutput>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithAlignedParaboloid(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AlignedParaboloidType& a_paraboloid, const double a_inv_volume_scale,
    SurfaceOutputType* a_surface = nullptr);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class AligneParaboloidType,
          class SurfaceOutputType = NoSurfaceOutput>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
formParaboloidIntersectionBases(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AligneParaboloidType& a_aligned_paraboloid,
    const UnsignedIndex_t a_nudge_iter, SurfaceOutputType* a_surface = nullptr);

template <class ScalarType>
inline NormalBase<ScalarType> computeTangentVectorAtPoint(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const NormalBase<ScalarType>& a_plane_normal,
    const PtBase<ScalarType>& a_pt);

template <class ScalarType, class PtWithGradientType>
inline PtWithGradientType computeTangentVectorAndGradientAtPoint(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtWithGradientType& a_plane_normal, const PtWithGradientType& a_pt);

template <class ScalarType>
inline NormalBase<ScalarType> computeAndCorrectTangentVectorAtPt(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const NormalBase<ScalarType>& a_plane_normal,
    const PtBase<ScalarType>& a_origin_pt, const PtBase<ScalarType>& a_end_pt,
    const NormalBase<ScalarType>& a_end_tangent,
    const PtBase<ScalarType>& a_intersection_pt);

template <class ScalarType, class PtTypeWithGradient>
inline PtTypeWithGradient computeAndCorrectTangentVectorAndGradientAtPt(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtTypeWithGradient& a_plane_normal,
    const PtTypeWithGradient& a_origin_pt, const PtTypeWithGradient& a_end_pt,
    const PtTypeWithGradient& a_end_tangent,
    const PtTypeWithGradient& a_intersection_pt);

template <class ReturnType, class ScalarType,
          class SurfaceOutputType = NoSurfaceOutput, class PtType>
ReturnType computeType3ContributionWithSplit(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const NormalBase<ScalarType>& a_plane_normal, const PtType& a_pt_ref,
    const PtType& a_pt_0, const PtType& a_pt_1,
    const NormalBase<ScalarType>& a_tangent_0,
    const NormalBase<ScalarType>& a_tangent_1, bool* a_requires_nudge,
    UnsignedIndex_t* a_split_counter, SurfaceOutputType* a_surface = nullptr);

template <class ReturnType, class ScalarType, class SurfaceOutputType,
          class PtType>
ReturnType computeType3ContributionWithGradientWithSplit(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtType& a_plane_normal, const PtType& a_pt_ref, const PtType& a_pt_0,
    const PtType& a_pt_1, const PtType& a_tangent_0, const PtType& a_tangent_1,
    SurfaceOutputType* a_surface);

template <class ReturnType, class ScalarType, class HalfEdgeType, class PtType>
ReturnType computeUnclippedSegmentType1Contribution(
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid,
    const PtType& a_ref_pt, const HalfEdgeType a_entry_half_edge,
    HalfEdgeType& a_exit_half_edge, const bool skip_first);

template <class ReturnType, class ScalarType, class HalfEdgeType,
          class SurfaceOutputType = NoSurfaceOutput, class PtType>
ReturnType computeNewEdgeSegmentContribution(
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid,
    const PtType& a_ref_pt, const HalfEdgeType a_entry_half_edge,
    const HalfEdgeType a_exit_half_edge, const bool skip_first,
    const bool a_ignore_type3, bool* a_requires_nudge,
    SurfaceOutputType* a_surface = nullptr);

template <class PtType, class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType>
void placeSingleIntercept(const PtType& a_intersection_location,
                          HalfEdgeType* a_half_edge_with_intersection,
                          SegmentedHalfEdgePolytopeType* a_polytope,
                          HalfEdgePolytopeType* a_complete_polytope);

template <class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType, class VertexType>
void placeDoubleIntercept(
    const StackVector<VertexType, 2>& a_intersection_location,
    HalfEdgeType* a_half_edge_with_intersection,
    SegmentedHalfEdgePolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope);

template <class PtType, class ScalarType>
void checkAndFindIntercepts(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid, const PtType& a_pt_0,
    const PtType& a_pt_1, StackVector<PtType, 2>* a_intercepts,
    const ScalarType a_nudge_epsilon, const bool a_elliptic);

template <class VertexType>
bool vertexBelow(
    const VertexType& a_pt,
    const AlignedParaboloidBase<typename VertexType::value_type>& a_paraboloid);

template <class ScalarType>
inline bool isPtBeforeIntersectionWithEdge(
    const std::array<ScalarType, 2>& a_test_pt,
    const PtBase<ScalarType>& a_vertex_0, const PtBase<ScalarType>& a_vertex_1);

template <class ScalarType>
inline bool isPtBeforeIntersectionWithEdgeWithComponent(
    const PtBase<ScalarType>& a_test_pt, const PtBase<ScalarType>& a_vertex_0,
    const PtBase<ScalarType>& a_vertex_1, const UnsignedIndex_t a_index);

template <class ScalarType, class HalfEdgeType>
bool ellipseContainedInFace(
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid,
    const PlaneBase<ScalarType>& a_face_plane, HalfEdgeType* const a_half_edge);

template <class ReturnType, class ScalarType, class HalfEdgeType,
          class SurfaceOutputType = NoSurfaceOutput>
ReturnType orientAndApplyType3Correction(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    HalfEdgeType* a_start, HalfEdgeType* a_end, bool* a_requires_nudge,
    SurfaceOutputType* a_surface = nullptr);

template <class ReturnType, class ScalarType, class HalfEdgeType,
          class SurfaceOutputType>
ReturnType orientAndApplyType3CorrectionWithGradients(
    const AlignedParaboloid& a_paraboloid, HalfEdgeType* a_start,
    HalfEdgeType* a_end, SurfaceOutputType* a_surface, bool* a_requires_nudge);

template <class ScalarType, class HalfEdgeType, class PtType>
Normal determineNudgeDirection(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    HalfEdgeType* a_current_edge, const PtType& a_inter_pt);

template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, void>
resetPolyhedron(SegmentedHalfEdgePolyhedronType* a_polytope,
                HalfEdgePolytopeType* a_complete_polytope);

template <class ScalarType>
AlignedParaboloidBase<Quad_t> nudgeParaboloid(
    AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const UnsignedIndex_t a_nudge_iter);

template <class ScalarType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class SurfaceOutputType>
void nudgePolyhedron(SegmentedHalfEdgePolyhedronType* a_polytope,
                     HalfEdgePolytopeType* a_complete_polytope,
                     const UnsignedIndex_t a_nudge_iter,
                     SurfaceOutputType* a_surface);

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class AligneParaboloidType,
          class SurfaceOutputType>
ReturnType reformParaboloidIntersectionBases(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AligneParaboloidType& a_aligned_paraboloid,
    const UnsignedIndex_t a_nudge_iter, SurfaceOutputType* a_surface);

template <class C>
struct has_embedded_gradient : std::false_type {};

template <class C>
struct has_embedded_gradient<const C> : has_embedded_gradient<C> {};

template <class GradientType>
struct has_embedded_gradient<PtWithGradient<GradientType>> : std::true_type {};

template <class GradientType>
struct has_embedded_gradient<VolumeWithGradient<GradientType>>
    : std::true_type {};

template <class GradientType>
struct has_embedded_gradient<VolumeMomentsWithGradient<GradientType>>
    : std::true_type {};

}  // namespace IRL

#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection.tpp"

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_H_
