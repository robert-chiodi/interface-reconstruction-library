// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_TPP_
#define IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_TPP_

#include <float.h>
#include <cassert>
#include <cmath>

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/generic_cutting/paraboloid_intersection/moment_contributions.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/helpers/mymath.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/paraboloid_reconstruction/rational_bezier_arc.h"

namespace IRL {

inline Normal computeTangentVectorAtPoint(const AlignedParaboloid& a_paraboloid,
                                          const Normal& a_plane_normal,
                                          const Pt& a_pt) {
  const Normal surface_normal = getParaboloidSurfaceNormal(a_paraboloid, a_pt);
  Normal tangent_at_pt = crossProduct(a_plane_normal, surface_normal);
  tangent_at_pt.normalize();
  return tangent_at_pt;
}

template <class PtWithGradientType>
inline PtWithGradientType computeTangentVectorAndGradientAtPoint(
    const AlignedParaboloid& a_paraboloid,
    const PtWithGradientType& a_plane_normal, const PtWithGradientType& a_pt) {
  using gradient_type = typename PtWithGradientType::gradient_type;
  const auto surface_normal_withgrad =
      getParaboloidSurfaceNormalWithGradient<PtWithGradientType>(a_paraboloid,
                                                                 a_pt);
  const Normal surface_normal = Normal::fromPt(surface_normal_withgrad.getPt());
  const auto& surface_normal_grad = surface_normal_withgrad.getData();
  const auto& a_plane_normal_grad = a_plane_normal.getData();
  const double tangent_at_pt_x = a_plane_normal[1] * surface_normal[2] -
                                 a_plane_normal[2] * surface_normal[1];
  const double tangent_at_pt_y = a_plane_normal[2] * surface_normal[0] -
                                 a_plane_normal[0] * surface_normal[2];
  const double tangent_at_pt_z = a_plane_normal[0] * surface_normal[1] -
                                 a_plane_normal[1] * surface_normal[0];
  const auto tangent_at_pt_gradx = a_plane_normal_grad[1] * surface_normal[2] +
                                   a_plane_normal[1] * surface_normal_grad[2] -
                                   a_plane_normal_grad[2] * surface_normal[1] -
                                   a_plane_normal[2] * surface_normal_grad[1];
  const auto tangent_at_pt_grady = a_plane_normal_grad[2] * surface_normal[0] +
                                   a_plane_normal[2] * surface_normal_grad[0] -
                                   a_plane_normal_grad[0] * surface_normal[2] -
                                   a_plane_normal[0] * surface_normal_grad[2];
  const auto tangent_at_pt_gradz = a_plane_normal_grad[0] * surface_normal[1] +
                                   a_plane_normal[0] * surface_normal_grad[1] -
                                   a_plane_normal_grad[1] * surface_normal[0] -
                                   a_plane_normal[1] * surface_normal_grad[0];
  const double norm_tangent = std::sqrt(tangent_at_pt_x * tangent_at_pt_x +
                                        tangent_at_pt_y * tangent_at_pt_y +
                                        tangent_at_pt_z * tangent_at_pt_z);
  const auto norm_tangent_grad = (tangent_at_pt_gradx * tangent_at_pt_x +
                                  tangent_at_pt_grady * tangent_at_pt_y +
                                  tangent_at_pt_gradz * tangent_at_pt_z) /
                                 safelyEpsilon(norm_tangent);
  assert(norm_tangent > DBL_EPSILON);
  const auto tangent_at_pt = Pt(tangent_at_pt_x / safelyEpsilon(norm_tangent),
                                tangent_at_pt_y / safelyEpsilon(norm_tangent),
                                tangent_at_pt_z / safelyEpsilon(norm_tangent));
  auto tangent_at_pt_withgrad = PtWithGradientType(tangent_at_pt);
  auto& tangent_at_pt_grad = tangent_at_pt_withgrad.getData();
  tangent_at_pt_grad[0] = (tangent_at_pt_gradx * norm_tangent -
                           tangent_at_pt_x * norm_tangent_grad) /
                          safelyEpsilon(norm_tangent * norm_tangent);
  tangent_at_pt_grad[1] = (tangent_at_pt_grady * norm_tangent -
                           tangent_at_pt_y * norm_tangent_grad) /
                          safelyEpsilon(norm_tangent * norm_tangent);
  tangent_at_pt_grad[2] = (tangent_at_pt_gradz * norm_tangent -
                           tangent_at_pt_z * norm_tangent_grad) /
                          safelyEpsilon(norm_tangent * norm_tangent);
  return tangent_at_pt_withgrad;
}

inline Normal computeAndCorrectTangentVectorAtPt(
    const AlignedParaboloid& a_paraboloid, const Normal& a_plane_normal,
    const Pt& a_origin_pt, const Pt& a_end_pt, const Normal& a_end_tangent,
    const Pt& a_intersection_pt) {
  Normal tangent = computeTangentVectorAtPoint(a_paraboloid, a_plane_normal,
                                               a_intersection_pt);
  const Normal correct_sign =
      crossProduct(a_plane_normal, a_end_pt - a_origin_pt);
  if ((a_end_tangent * correct_sign > 0.0) !=
      (tangent * Normal(crossProduct(a_plane_normal,
                                     a_intersection_pt - a_origin_pt)) >
       0.0)) {
    tangent = -tangent;
  }
  return tangent;
}

template <class PtTypeWithGradient>
inline PtTypeWithGradient computeAndCorrectTangentVectorAndGradientAtPt(
    const AlignedParaboloid& a_paraboloid,
    const PtTypeWithGradient& a_plane_normal,
    const PtTypeWithGradient& a_origin_pt, const PtTypeWithGradient& a_end_pt,
    const PtTypeWithGradient& a_end_tangent,
    const PtTypeWithGradient& a_intersection_pt) {
  PtTypeWithGradient tangent_withgrad = computeTangentVectorAndGradientAtPoint(
      a_paraboloid, a_plane_normal, a_intersection_pt);
  const Normal tangent = Normal::fromPt(tangent_withgrad.getPt());
  const Normal plane_normal = Normal::fromPt(a_plane_normal.getPt());
  const Normal edge_vector1 =
      Normal::fromPt(a_end_pt.getPt() - a_origin_pt.getPt());
  const Normal edge_vector2 =
      Normal::fromPt(a_intersection_pt.getPt() - a_origin_pt.getPt());
  const Normal end_tangent = Normal::fromPt(a_end_tangent.getPt());
  const Normal correct_sign = crossProduct(plane_normal, edge_vector1);
  if ((end_tangent * correct_sign > 0.0) !=
      (tangent * Normal(crossProduct(plane_normal, edge_vector2)) > 0.0)) {
    tangent_withgrad = -tangent_withgrad;
  }
  return tangent_withgrad;
}

template <class ReturnType, class SurfaceOutputType, class PtType>
ReturnType computeType3ContributionWithSplit(
    const AlignedParaboloid& a_paraboloid, const Normal& a_plane_normal,
    const PtType& a_pt_ref, const PtType& a_pt_0, const PtType& a_pt_1,
    const Normal& a_tangent_0, const Normal& a_tangent_1,
    SurfaceOutputType* a_surface) {
  const Pt& pt_ref = a_pt_ref.getPt();
  const Pt& pt_0 = a_pt_0.getPt();
  const Pt& pt_1 = a_pt_1.getPt();
  if (squaredMagnitude(pt_1 - pt_0) < 100.0 * DBL_EPSILON * DBL_EPSILON) {
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      a_surface->addArc(
          RationalBezierArc(pt_0, 0.5 * (pt_0 + pt_1), pt_1, 0.0));
    }
    return ReturnType::fromScalarConstant(0.0);
  } else {
    const Normal edge_vector = pt_1 - pt_0;
    bool split = ((a_tangent_0 * edge_vector) < 0.0 ||
                  (a_tangent_1 * edge_vector) > 0.0);
    if (split) {
      const Pt average_pt = 0.5 * (pt_0 + pt_1);
      auto average_tangent = Normal(0.5 * (a_tangent_0 + a_tangent_1));
      if (squaredMagnitude(average_tangent) <
          100.0 * DBL_EPSILON * DBL_EPSILON) {
        average_tangent = Normal(0.25 * a_tangent_0 + 0.75 * a_tangent_1);
      }
      average_tangent.normalize();
      Pt projected_pt = projectPtAlongHalfLineOntoParaboloid(
          a_paraboloid, average_tangent, average_pt);
      const Normal tangent_projected_pt = computeAndCorrectTangentVectorAtPt(
          a_paraboloid, a_plane_normal, pt_0, pt_1, a_tangent_1, projected_pt);
      // We need to store this vertex so that its address remains unique over
      // time
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        Pt* new_point = new Pt(projected_pt);
        a_surface->addPt(new_point);
        return computeType3ContributionWithSplit<ReturnType>(
                   a_paraboloid, a_plane_normal, a_pt_ref, a_pt_0, *new_point,
                   a_tangent_0, tangent_projected_pt, a_surface) +
               computeType3ContributionWithSplit<ReturnType>(
                   a_paraboloid, a_plane_normal, a_pt_ref, *new_point, a_pt_1,
                   -tangent_projected_pt, a_tangent_1, a_surface);
      } else {
        return computeType3ContributionWithSplit<ReturnType>(
                   a_paraboloid, a_plane_normal, a_pt_ref, a_pt_0,
                   PtType(projected_pt), a_tangent_0, tangent_projected_pt,
                   a_surface) +
               computeType3ContributionWithSplit<ReturnType>(
                   a_paraboloid, a_plane_normal, a_pt_ref, PtType(projected_pt),
                   a_pt_1, -tangent_projected_pt, a_tangent_1, a_surface);
      }
    } else {
      const auto arc = RationalBezierArc(pt_0, a_tangent_0, pt_1, a_tangent_1,
                                         a_plane_normal, a_paraboloid);
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        a_surface->addArc(arc);
      }

      auto moments = computeType3Contribution<ReturnType>(a_paraboloid, arc);
      if (!(&a_pt_ref == &a_pt_0 || &a_pt_ref == &a_pt_1)) {
        moments += computeTriangleCorrection<ReturnType>(a_paraboloid, pt_0,
                                                         pt_1, pt_ref);
      }
      return moments;
    }
  }
}

template <class ReturnType, class SurfaceOutputType, class PtType>
ReturnType computeType3ContributionWithGradientWithSplit(
    const AlignedParaboloid& a_paraboloid, const PtType& a_plane_normal,
    const PtType& a_pt_ref, const PtType& a_pt_0, const PtType& a_pt_1,
    const PtType& a_tangent_0, const PtType& a_tangent_1,
    SurfaceOutputType* a_surface) {
  using gradient_type = typename PtType::gradient_type;
  const Pt& pt_ref = a_pt_ref.getPt();
  const Pt& pt_0 = a_pt_0.getPt();
  const Pt& pt_1 = a_pt_1.getPt();
  const Normal tgt_0 = Normal::fromPt(a_tangent_0.getPt());
  const Normal tgt_1 = Normal::fromPt(a_tangent_1.getPt());
  if (squaredMagnitude(pt_1 - pt_0) < 100.0 * DBL_EPSILON * DBL_EPSILON) {
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      a_surface->addArc(
          RationalBezierArc(pt_0, 0.5 * (pt_0 + pt_1), pt_1, 0.0));
    }
    return ReturnType::fromScalarConstant(0.0);
  } else {
    const Normal edge_vector = pt_1 - pt_0;
    bool split = ((tgt_0 * edge_vector) < 0.0 || (tgt_1 * edge_vector) > 0.0);
    if (split) {
      const auto average_pt_withgrad = 0.5 * (a_pt_0 + a_pt_1);
      auto avg_tgt_withgrad = 0.5 * (a_tangent_0 + a_tangent_1);
      if (squaredMagnitude(avg_tgt_withgrad.getPt()) <
          100.0 * DBL_EPSILON * DBL_EPSILON) {
        avg_tgt_withgrad = 0.25 * a_tangent_0 + 0.75 * a_tangent_1;
      }
      const Pt avg_tgt = avg_tgt_withgrad.getPt();
      const auto avg_tgt_grad = avg_tgt_withgrad.getData();
      const double norm_avg_tgt = magnitude(avg_tgt);
      const gradient_type norm_avg_tgt_grad =
          (avg_tgt_grad[0] * avg_tgt[0] + avg_tgt_grad[1] * avg_tgt[1] +
           avg_tgt_grad[2] * avg_tgt[2]) /
          safelyEpsilon(norm_avg_tgt);
      auto avg_normalized_tgt_withgrad =
          PtType(avg_tgt / safelyEpsilon(norm_avg_tgt));
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        avg_normalized_tgt_withgrad.getData()[d] =
            (avg_tgt_grad[d] * norm_avg_tgt - avg_tgt[d] * norm_avg_tgt_grad) /
            safelyEpsilon(norm_avg_tgt * norm_avg_tgt);
      }
      PtType projected_pt_withgrad =
          projectPtAlongHalfLineOntoParaboloidWithGradient<PtType>(
              a_paraboloid, avg_tgt_withgrad, average_pt_withgrad);
      const auto tangent_projected_pt_withgrad =
          computeAndCorrectTangentVectorAndGradientAtPt<PtType>(
              a_paraboloid, a_plane_normal, a_pt_0, a_pt_1, a_tangent_1,
              projected_pt_withgrad);
      // We need to store this vertex so that its address remains unique over
      // time
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        Pt* new_point = new Pt(projected_pt_withgrad.getPt());
        a_surface->addPt(new_point);
        auto new_proj_pt_withgrad = PtType(*new_point);
        new_proj_pt_withgrad.getData() = projected_pt_withgrad.getData();
        return computeType3ContributionWithGradientWithSplit<ReturnType>(
                   a_paraboloid, a_plane_normal, a_pt_ref, a_pt_0,
                   new_proj_pt_withgrad, a_tangent_0,
                   tangent_projected_pt_withgrad, a_surface) +
               computeType3ContributionWithGradientWithSplit<ReturnType>(
                   a_paraboloid, a_plane_normal, a_pt_ref, new_proj_pt_withgrad,
                   a_pt_1, -tangent_projected_pt_withgrad, a_tangent_1,
                   a_surface);
      } else {
        return computeType3ContributionWithGradientWithSplit<ReturnType>(
                   a_paraboloid, a_plane_normal, a_pt_ref, a_pt_0,
                   projected_pt_withgrad, a_tangent_0,
                   tangent_projected_pt_withgrad, a_surface) +
               computeType3ContributionWithGradientWithSplit<ReturnType>(
                   a_paraboloid, a_plane_normal, a_pt_ref,
                   projected_pt_withgrad, a_pt_1,
                   -tangent_projected_pt_withgrad, a_tangent_1, a_surface);
      }
    } else {
      const auto arc_with_gradient = RationalBezierArcWithGradient<PtType>(
          a_pt_0, a_tangent_0, a_pt_1, a_tangent_1, a_plane_normal,
          a_paraboloid);
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        const auto arc = RationalBezierArc(
            a_pt_0.getPt(), arc_with_gradient.control_point().getPt(),
            a_pt_1.getPt(), arc_with_gradient.weight());
        a_surface->addArc(arc);
      }
      auto volume_with_gradient =
          computeType3Contribution<ReturnType>(a_paraboloid, arc_with_gradient);
      if (!(&a_pt_ref == &a_pt_0 || &a_pt_ref == &a_pt_1)) {
        volume_with_gradient += computeTriangleCorrection<ReturnType>(
            a_paraboloid, a_pt_0, a_pt_1, a_pt_ref);
      }
      return volume_with_gradient;
    }
  }
}

// Starts from an entry, returns the exit that is reached.
template <class ReturnType, class HalfEdgeType, class PtType>
ReturnType computeUnclippedSegmentType1Contribution(
    const AlignedParaboloid& a_aligned_paraboloid, const PtType& a_ref_pt,
    const HalfEdgeType a_entry_half_edge, HalfEdgeType& a_exit_half_edge,
    const bool skip_first) {
  ReturnType full_moments = ReturnType::fromScalarConstant(0.0);

  assert(a_entry_half_edge->getPreviousVertex()->isClipped() ||
         (a_entry_half_edge->getPreviousVertex()->doesNotNeedToSeek() &&
          a_entry_half_edge->getNextHalfEdge()->getVertex()->isNotClipped()));

  auto current_half_edge = a_entry_half_edge->getNextHalfEdge();
  auto vertex = current_half_edge->getVertex();
  if (skip_first && vertex->needsToSeek()) {
    current_half_edge = current_half_edge->getNextHalfEdge();
  }
  auto prev_pt = current_half_edge->getPreviousVertex()->getLocation();

  if (!(skip_first && vertex->doesNotNeedToSeek())) {
    while (true) {
      vertex = current_half_edge->getVertex();
      const auto& curr_pt = vertex->getLocation();
      full_moments +=
          computeType1Contribution<ReturnType>(a_ref_pt, prev_pt, curr_pt);
      if (vertex->needsToSeek()) {
        prev_pt = curr_pt;
        current_half_edge = current_half_edge->getNextHalfEdge();
      } else {
        assert(!(current_half_edge->getPreviousVertex()->isClipped() ||
                 (current_half_edge->getPreviousVertex()->doesNotNeedToSeek() &&
                  current_half_edge->getNextHalfEdge()
                      ->getVertex()
                      ->isNotClipped())));
        break;
      }
    }
  }
  a_exit_half_edge = current_half_edge;
  return full_moments;
}

template <class ReturnType, class HalfEdgeType, class SurfaceOutputType,
          class PtType>
ReturnType computeNewEdgeSegmentContribution(
    const AlignedParaboloid& a_aligned_paraboloid, const PtType& a_ref_pt,
    const HalfEdgeType a_entry_half_edge, const HalfEdgeType a_exit_half_edge,
    const bool skip_first, SurfaceOutputType* a_surface) {
  ReturnType full_moments = ReturnType::fromScalarConstant(0.0);
  // Handle new edge on exit->entry
  if (!skip_first) {
    full_moments += computeType1Contribution<ReturnType>(
        a_ref_pt, a_exit_half_edge->getVertex()->getLocation(),
        a_entry_half_edge->getVertex()->getLocation());
  }

  full_moments += computeType2Contribution<ReturnType>(
      a_aligned_paraboloid, a_exit_half_edge->getVertex()->getLocation(),
      a_entry_half_edge->getVertex()->getLocation());

  if constexpr (has_embedded_gradient<ReturnType>::value) {
    full_moments += orientAndApplyType3CorrectionWithGradients<ReturnType>(
        a_aligned_paraboloid, a_exit_half_edge, a_entry_half_edge, a_surface);
  } else {
    full_moments += orientAndApplyType3Correction<ReturnType>(
        a_aligned_paraboloid, a_exit_half_edge, a_entry_half_edge, a_surface);
  }
  return full_moments;
}

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithParaboloid(SegmentedHalfEdgePolyhedronType* a_polytope,
                                  HalfEdgePolytopeType* a_complete_polytope,
                                  const Paraboloid& a_paraboloid) {
  ReturnType moments;
  if (a_paraboloid.isAlwaysAbove()) {
    if constexpr (has_paraboloid_surface<ReturnType>::value) {
      moments.getMoments() =
          ReturnType::moment_type::calculateMoments(a_polytope);
      moments.getSurface().setParaboloid(a_paraboloid);
    } else {
      moments = ReturnType::calculateMoments(a_polytope);
    }
    return moments;
  } else if (a_paraboloid.isAlwaysBelow()) {
    if constexpr (has_paraboloid_surface<ReturnType>::value) {
      moments.getMoments() = 0.0;
      moments.getSurface().setParaboloid(a_paraboloid);
    } else {
      moments = 0.0;
    }
    return moments;
  }

  // Move into reconstruction reference frame
  const UnsignedIndex_t original_number_of_vertices =
      a_polytope->getNumberOfVertices();
  const auto& datum = a_paraboloid.getDatum();
  const auto& ref_frame = a_paraboloid.getReferenceFrame();
  assert(ref_frame.isOrthonormalBasis());

  for (UnsignedIndex_t v = 0; v < original_number_of_vertices; ++v) {
    const Pt original_pt =
        a_polytope->getVertex(v)->getLocation().getPt() - datum;
    typename SegmentedHalfEdgePolyhedronType::pt_type projected_location;
    auto& pt = projected_location.getPt();
    for (UnsignedIndex_t n = 0; n < 3; ++n) {
      pt[n] = ref_frame[n] * original_pt;
    }
    a_polytope->getVertex(v)->setLocation(projected_location);
  }

  // Recalculate in face plane information
  for (auto& face : (*a_polytope)) {
    auto normal = Normal(0.0, 0.0, 0.0);
    const auto starting_half_edge = face->getStartingHalfEdge();
    auto current_half_edge = starting_half_edge;
    auto next_half_edge = starting_half_edge->getNextHalfEdge();
    const auto& start_location =
        starting_half_edge->getPreviousVertex()->getLocation().getPt();
    do {
      normal += crossProduct(
          current_half_edge->getVertex()->getLocation().getPt() -
              start_location,
          next_half_edge->getVertex()->getLocation().getPt() - start_location);
      current_half_edge = next_half_edge;
      next_half_edge = next_half_edge->getNextHalfEdge();
    } while (next_half_edge != starting_half_edge);
    normal.normalize();
    face->setPlane(Plane(normal, normal * start_location));
  }

  // Compute intersection
  if constexpr (has_paraboloid_surface<ReturnType>::value) {
    moments.getSurface().setParaboloid(a_paraboloid);
    moments.getMoments() =
        intersectPolyhedronWithParaboloid<typename ReturnType::moment_type>(
            a_polytope, a_complete_polytope,
            a_paraboloid.getAlignedParaboloid(), &moments.getSurface());
  } else {
    NoSurfaceOutput* surf = nullptr;
    moments = intersectPolyhedronWithParaboloid<ReturnType>(
        a_polytope, a_complete_polytope, a_paraboloid.getAlignedParaboloid(),
        surf);
  }

  // Rotate base polyhedron back
  const UnsignedIndex_t number_of_vertices = a_polytope->getNumberOfVertices();
  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    const Pt& original_pt = a_polytope->getVertex(v)->getLocation().getPt();
    typename SegmentedHalfEdgePolyhedronType::pt_type projected_location;
    auto& pt = projected_location.getPt();
    pt = Pt(0.0, 0.0, 0.0);
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      for (UnsignedIndex_t n = 0; n < 3; ++n) {
        pt[n] += ref_frame[d][n] * original_pt[d];
      }
    }
    pt += datum;
    a_polytope->getVertex(v)->setLocation(projected_location);
  }

  // Move first moment back to original frame of reference
  if constexpr (has_paraboloid_surface<ReturnType>::value) {
    if constexpr (std::is_same<typename ReturnType::moment_type,
                               VolumeMoments>::value) {
      auto pt = Pt(0.0, 0.0, 0.0);
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          pt[n] += ref_frame[d][n] * moments.centroid()[d];
        }
      }
      pt += datum;
      moments.centroid() = pt;
    }
  }
  if constexpr (!has_paraboloid_surface<ReturnType>::value &&
                std::is_same<ReturnType, VolumeMoments>::value) {
    auto pt = Pt(0.0, 0.0, 0.0);
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      for (UnsignedIndex_t n = 0; n < 3; ++n) {
        pt[n] += ref_frame[d][n] * moments.centroid()[d];
      }
    }
    pt += datum;
    moments.centroid() = pt;
  }

  return moments;
}

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class SurfaceOutputType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithParaboloid(SegmentedHalfEdgePolyhedronType* a_polytope,
                                  HalfEdgePolytopeType* a_complete_polytope,
                                  const AlignedParaboloid& a_paraboloid,
                                  SurfaceOutputType* a_surface) {
  // Assumed a_polytope is already rotated to be in same coordinate system
  // as a_paraboloid.

  // Below function computes the entire integration
  const auto moments = formParaboloidIntersectionBases<ReturnType>(
      a_polytope, a_complete_polytope, a_paraboloid, 0, a_surface);

  return moments;
}

template <class PtType, class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType>
void placeSingleIntercept(const PtType& a_intersection_location,
                          HalfEdgeType* a_half_edge_with_intersection,
                          SegmentedHalfEdgePolytopeType* a_polytope,
                          HalfEdgePolytopeType* a_complete_polytope) {
  auto first_intersection_vertex = a_complete_polytope->getNewVertex(
      typename SegmentedHalfEdgePolytopeType::vertex_type(
          a_intersection_location));
  first_intersection_vertex->setToSeek();
  a_polytope->addVertex(first_intersection_vertex);
  HalfEdgeType* new_half_edge = separateIntersectedHalfEdge(
      first_intersection_vertex, a_half_edge_with_intersection, a_polytope,
      a_complete_polytope);
  createOppositeHalfEdgeFromIntersection(a_half_edge_with_intersection,
                                         new_half_edge, a_complete_polytope);
}

template <class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType, class VertexType>
void placeDoubleIntercept(const StackVector<std::pair<VertexType, double>, 2>&
                              a_intersection_location,
                          HalfEdgeType* a_half_edge_with_intersection,
                          SegmentedHalfEdgePolytopeType* a_polytope,
                          HalfEdgePolytopeType* a_complete_polytope) {
  assert(a_intersection_location.size() == 2);

  // Need to place furthest vertex first so that a_half_edge_with_intersection
  // remains attached to same vertex it started on. Also want to keep property
  // that the half edge the vertex stores has the previousVertex also
  // unclipped. This means for the first (furthest) vertex, need to reference
  // the opposite half edge of the current one. This enables more efficient
  // face truncation in the next phase of the algorithm, since new vertices
  // will always have a half edge going from unclipped->new->clipped.

  placeSingleIntercept(a_intersection_location[0].first,
                       a_half_edge_with_intersection, a_polytope,
                       a_complete_polytope);
  a_half_edge_with_intersection->getPreviousVertex()->setHalfEdge(
      a_half_edge_with_intersection->getOppositeHalfEdge());
  placeSingleIntercept(a_intersection_location[1].first,
                       a_half_edge_with_intersection, a_polytope,
                       a_complete_polytope);
}

template <class PtType>
void checkAndFindIntercepts(
    const AlignedParaboloid& a_paraboloid, const PtType& a_pt_0,
    const PtType& a_pt_1,
    StackVector<std::pair<PtType, double>, 2>* a_intercepts,
    const double a_nudge_epsilon);

template <>
inline enable_if_t<!has_embedded_gradient<Pt>::value, void>
checkAndFindIntercepts(const AlignedParaboloid& a_paraboloid, const Pt& a_pt_0,
                       const Pt& a_pt_1,
                       StackVector<std::pair<Pt, double>, 2>* a_intercepts,
                       const double a_nudge_epsilon) {
  const auto& pt_0 = a_pt_0.getPt();
  const auto& pt_1 = a_pt_1.getPt();
  const auto pt_diff = pt_1 - pt_0;
  const double a = a_paraboloid.a() * pt_diff[0] * pt_diff[0] +
                   a_paraboloid.b() * pt_diff[1] * pt_diff[1];
  const double b = 2.0 * (a_paraboloid.a() * pt_diff[0] * pt_0[0] +
                          a_paraboloid.b() * pt_diff[1] * pt_0[1]) +
                   pt_diff[2];
  const double c = a_paraboloid.a() * pt_0[0] * pt_0[0] +
                   a_paraboloid.b() * pt_0[1] * pt_0[1] + pt_0[2];

  a_intercepts->resize(0);
  const auto solutions = solveQuadratic(a, b, c);
  for (const auto& solution : solutions) {
    if (solution >= -a_nudge_epsilon && solution <= 1.0 + a_nudge_epsilon) {
      a_intercepts->push_back(
          std::pair<Pt, double>(Pt(pt_0 + solution * pt_diff), solution));
    }
  }
}

template <class PtType>
inline enable_if_t<has_embedded_gradient<PtType>::value, void>
checkAndFindIntercepts(
    const AlignedParaboloid& a_paraboloid,
    const PtWithGradient<typename PtType::gradient_type>& a_pt_0,
    const PtWithGradient<typename PtType::gradient_type>& a_pt_1,
    StackVector<
        std::pair<PtWithGradient<typename PtType::gradient_type>, double>, 2>*
        a_intercepts,
    const double a_nudge_epsilon) {
  using gradient_type = typename PtType::gradient_type;
  const auto pt_diff_with_grad = a_pt_1 - a_pt_0;
  const auto& pt_0 = a_pt_0.getPt();
  const auto& pt_1 = a_pt_1.getPt();
  const auto& pt_diff = pt_diff_with_grad.getPt();
  const auto& pt_0_grad = a_pt_0.getData();
  const auto& pt_diff_grad = pt_diff_with_grad.getData();
  const double A = a_paraboloid.a(), B = a_paraboloid.b();
  auto A_grad = gradient_type(0.0), B_grad = gradient_type(0.0);
  A_grad.setGradA(1.0);
  B_grad.setGradB(1.0);
  const double a = A * pt_diff[0] * pt_diff[0] + B * pt_diff[1] * pt_diff[1];
  const double b =
      2.0 * (A * pt_diff[0] * pt_0[0] + B * pt_diff[1] * pt_0[1]) + pt_diff[2];
  const double c = A * pt_0[0] * pt_0[0] + B * pt_0[1] * pt_0[1] + pt_0[2];
  const auto a_grad = A_grad * pt_diff[0] * pt_diff[0] +
                      B_grad * pt_diff[1] * pt_diff[1] +
                      A * (2.0 * pt_diff_grad[0] * pt_diff[0]) +
                      B * (2.0 * pt_diff_grad[1] * pt_diff[1]);
  const auto b_grad =
      2.0 * (A_grad * pt_diff[0] * pt_0[0] + B_grad * pt_diff[1] * pt_0[1] +
             A * (pt_diff_grad[0] * pt_0[0] + pt_diff[0] * pt_0_grad[0]) +
             B * (pt_diff_grad[1] * pt_0[1] + pt_diff[1] * pt_0_grad[1])) +
      pt_diff_grad[2];
  const auto c_grad = A_grad * pt_0[0] * pt_0[0] + B_grad * pt_0[1] * pt_0[1] +
                      A * 2.0 * pt_0_grad[0] * pt_0[0] +
                      B * 2.0 * pt_0_grad[1] * pt_0[1] + pt_0_grad[2];
  a_intercepts->resize(0);
  const auto solutions = solveQuadraticWithGradient<gradient_type>(
      a, b, c, a_grad, b_grad, c_grad);
  for (const auto& solution : solutions) {
    const auto sol = solution.first;
    const auto sol_grad = solution.second;
    if (sol >= -a_nudge_epsilon && sol <= 1.0 + a_nudge_epsilon) {
      auto intersection = PtType(Pt(pt_0 + sol * pt_diff));
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        intersection.getData()[d] =
            pt_0_grad[d] + sol_grad * pt_diff[d] + sol * pt_diff_grad[d];
      }
      a_intercepts->push_back(std::pair<PtType, double>({intersection, sol}));
    }
  }
}

template <class VertexType>
bool vertexBelow(const VertexType& a_pt,
                 const AlignedParaboloid& a_paraboloid) {
  const auto& pt = a_pt.getPt();
  return pt[2] <
         -(a_paraboloid.a() * pt[0] * pt[0] + a_paraboloid.b() * pt[1] * pt[1]);
}

// This algorithm is based on the method shared from the website below, and
// from several stackoverflow posts citing that website.
// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
inline bool isPtBeforeIntersectionWithEdge(
    const std::array<double, 2>& a_test_pt, const Pt& a_vertex_0,
    const Pt& a_vertex_1) {
  if ((a_test_pt[1] > a_vertex_0[1]) == (a_test_pt[1] > a_vertex_1[1])) {
    return false;  // Projected ray never intersects edge.
  }
  double location_of_intersection_along_ray =
      (a_vertex_0[0] - a_vertex_1[0]) * (a_test_pt[1] - a_vertex_1[1]) /
          (a_vertex_0[1] - a_vertex_1[1]) +
      a_vertex_1[0];
  // Intersection was to the right if location_of_intersection_along_ray is
  // greater.
  return a_test_pt[0] < location_of_intersection_along_ray;
}

// This algorithm is based on the method shared from the website below, and
// from several stackoverflow posts citing that website.
// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
inline bool isPtBeforeIntersectionWithEdgeWithComponent(
    const Pt& a_test_pt, const Pt& a_vertex_0, const Pt& a_vertex_1,
    const UnsignedIndex_t a_index) {
  const UnsignedIndex_t id0 = a_index;
  const UnsignedIndex_t id1 = (a_index + 1) % 3;
  if ((a_test_pt[id1] > a_vertex_0[id1]) ==
      (a_test_pt[id1] > a_vertex_1[id1])) {
    return false;  // Projected ray never intersects edge.
  }
  double location_of_intersection_along_ray =
      (a_vertex_0[id0] - a_vertex_1[id0]) * (a_test_pt[id1] - a_vertex_1[id1]) /
          (a_vertex_0[id1] - a_vertex_1[id1]) +
      a_vertex_1[id0];
  // Intersection was to the right if location_of_intersection_along_ray is
  // greater.
  return a_test_pt[id0] < location_of_intersection_along_ray;
}

// If centroid is outside of polygon, or any vertex on face
// is inside the ellipse, then the ellipse is not contained by the face.
template <class HalfEdgeType>
bool ellipseContainedInFace(const AlignedParaboloid& a_aligned_paraboloid,
                            const Plane& a_face_plane,
                            HalfEdgeType* const a_half_edge) {
  const auto& face_normal = a_face_plane.normal();
  const std::array<double, 2> conic_center{
      {face_normal[0] / (2.0 * a_aligned_paraboloid.a() * face_normal[2]),
       face_normal[1] / (2.0 * a_aligned_paraboloid.b() * face_normal[2])}};
  const double delta_face = a_face_plane.distance() / face_normal[2];
  const double gamma_face =
      a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
      a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] - delta_face;
  if (a_aligned_paraboloid.a() * gamma_face < 0.0) {
    return false;
  }

  // First we will check if centroid is in the bounding box
  // of the face polygon. Due to the fact we are assuming
  // the paraboloid was in Z direction, we assume the ellipse
  // lives on the x/y plane and project the polygon
  // down to it as well (essentially neglecting the z component).
  auto current_half_edge = a_half_edge;
  std::array<double, 2> xy_min{{DBL_MAX, DBL_MAX}};
  std::array<double, 2> xy_max{{-DBL_MAX, -DBL_MAX}};
  do {
    const Pt& location = current_half_edge->getVertex()->getLocation().getPt();
    for (UnsignedIndex_t d = 0; d < 2; ++d) {
      xy_min[d] = std::min(xy_min[d], location[d]);
      xy_max[d] = std::max(xy_max[d], location[d]);
    }
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != a_half_edge);
  if (conic_center[0] < xy_min[0] || conic_center[0] > xy_max[0] ||
      conic_center[1] < xy_min[1] || conic_center[1] > xy_max[1]) {
    return false;
  }

  // If not outside bounding box, need to use
  // a more refined test based on the number
  // of edges a ray emitted from the point crosses.
  // Even edges means the point is outside the face.
  current_half_edge = a_half_edge;
  bool pt_internal_to_polygon = false;
  do {
    const Pt& location_0 =
        current_half_edge->getPreviousVertex()->getLocation().getPt();
    const Pt& location_1 =
        current_half_edge->getVertex()->getLocation().getPt();
    if (isPtBeforeIntersectionWithEdge(conic_center, location_0, location_1)) {
      pt_internal_to_polygon = !pt_internal_to_polygon;
    }
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != a_half_edge);
  return pt_internal_to_polygon;
}

template <class ReturnType, class HalfEdgeType, class SurfaceOutputType>
ReturnType orientAndApplyType3Correction(const AlignedParaboloid& a_paraboloid,
                                         HalfEdgeType* a_start,
                                         HalfEdgeType* a_end,
                                         SurfaceOutputType* a_surface) {
  const auto& pt_0 = a_start->getVertex()->getLocation().getPt();
  const auto& pt_1 = a_end->getVertex()->getLocation().getPt();
  const auto edge_vector = Normal(pt_1 - pt_0);
  const auto& face_plane = a_end->getFace()->getPlane();
  const auto& face_normal = face_plane.normal();
  Normal tgt_0 = computeTangentVectorAtPoint(a_paraboloid, face_normal, pt_0);
  Normal tgt_1 = computeTangentVectorAtPoint(a_paraboloid, face_normal, pt_1);
  const bool elliptic_face =
      a_paraboloid.a() * a_paraboloid.b() > 0.0 && face_normal[2] != 0.0;

  if (!elliptic_face)  // The arc is from a hyperbola or parabola
  {
    const Normal n_cross_t0 = crossProduct(face_normal, tgt_0);
    if (std::fabs(n_cross_t0 * tgt_1) <
        100.0 * DBL_EPSILON)  // The tangents are (close to being) parallel
    {
      // We consider an arc with zero contribution
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        a_surface->addArc(
            RationalBezierArc(pt_1, 0.5 * (pt_0 + pt_1), pt_0, 0.0));
      }
      return ReturnType::fromScalarConstant(0.0);
    } else {  // The tangents are NOT parallel
      // Compute control point
      const double lambda_1 =
          -(n_cross_t0 * edge_vector) / (n_cross_t0 * tgt_1);
      const auto control_pt = Pt(pt_1 + lambda_1 * tgt_1);
      // Orient tangents to point towards control point
      tgt_0 = (tgt_0 * Normal(control_pt - pt_0)) < 0.0 ? -tgt_0 : tgt_0;
      tgt_1 = (tgt_1 * Normal(control_pt - pt_1)) < 0.0 ? -tgt_1 : tgt_1;
      // Construct Bezier arc from tangents
      const auto arc = RationalBezierArc(pt_1, tgt_1, pt_0, tgt_0, face_normal,
                                         a_paraboloid);
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        a_surface->addArc(arc);
      }
      return computeType3Contribution<ReturnType>(a_paraboloid, arc);
    }
  } else  // The arc is from an ellipse
  {
    // Compute edge vectors and check if tangent is parallel to edge
    Normal edge_00 = a_start->getVertex()->getLocation().getPt() -
                     a_start->getPreviousVertex()->getLocation().getPt();
    Normal edge_01 =
        a_start->getNextHalfEdge()->getVertex()->getLocation().getPt() -
        a_start->getVertex()->getLocation().getPt();
    Normal edge_0 = (squaredMagnitude(edge_00) > squaredMagnitude(edge_01))
                        ? edge_00
                        : edge_01;
    edge_0.normalize();
    Normal edge_10 = a_end->getVertex()->getLocation().getPt() -
                     a_end->getPreviousVertex()->getLocation().getPt();
    Normal edge_11 =
        a_end->getNextHalfEdge()->getVertex()->getLocation().getPt() -
        a_end->getVertex()->getLocation().getPt();
    Normal edge_1 = (squaredMagnitude(edge_10) > squaredMagnitude(edge_11))
                        ? edge_10
                        : edge_11;
    edge_1.normalize();
    bool tgt_0_parallel_edge_0 =
        std::abs(1.0 - std::fabs(tgt_0 * edge_0)) < 100.0 * DBL_EPSILON;
    bool tgt_1_parallel_edge_1 =
        std::abs(1.0 - std::fabs(tgt_1 * edge_1)) < 100.0 * DBL_EPSILON;
    if (!tgt_0_parallel_edge_0 &&
        !tgt_1_parallel_edge_1)  // We orient the tangent with the edge
                                 // normal
    {
      const Normal edge_normal_0 = crossProduct(face_normal, edge_0);
      const Normal edge_normal_1 = crossProduct(face_normal, edge_1);
      tgt_0 = (edge_normal_0 * tgt_0 < 0.0) ? -tgt_0 : tgt_0;
      tgt_1 = (edge_normal_1 * tgt_1 < 0.0) ? -tgt_1 : tgt_1;
    } else {
      // Compute ellipse center an orient tangents accordingly
      const Pt conic_center = conicCenter(face_plane, a_paraboloid);
      auto center_to_pt_0 = Normal(pt_0 - conic_center);
      center_to_pt_0.normalize();
      auto center_to_pt_1 = Normal(pt_1 - conic_center);
      center_to_pt_1.normalize();
      Normal dummy_tgt_0 = crossProduct(face_normal, center_to_pt_0);
      Normal dummy_tgt_1 = crossProduct(face_normal, center_to_pt_1);
      assert(std::abs(tgt_0 * dummy_tgt_0) > 10.0 * DBL_EPSILON);
      assert(std::abs(tgt_1 * dummy_tgt_1) > 10.0 * DBL_EPSILON);
      tgt_0 = (tgt_0 * dummy_tgt_0) < 0.0 ? -tgt_0 : tgt_0;
      tgt_1 = (tgt_1 * dummy_tgt_1) > 0.0 ? -tgt_1 : tgt_1;
      // At this point, the tangents form a valid arc (but they may be
      // oriented in the wrong direction)
      if (!tgt_0_parallel_edge_0) {
        const Normal edge_normal_0 = crossProduct(face_normal, edge_0);
        if (edge_normal_0 * tgt_0 < 0.0) {
          tgt_0 = -tgt_0;
          tgt_1 = -tgt_1;
        }
      } else if (!tgt_1_parallel_edge_1) {
        const Normal edge_normal_1 = crossProduct(face_normal, edge_1);
        if (edge_normal_1 * tgt_1 < 0.0) {
          tgt_0 = -tgt_0;
          tgt_1 = -tgt_1;
        }
      } else {
        // The arcs should have been sorted in this case!! This allows to
        // project the mid-point onto the arc, and check that it belongs to
        // the face
        const Pt avg_pt = 0.5 * (pt_0 + pt_1);
        auto avg_tgt = Normal(0.5 * (tgt_0 + tgt_1));
        if (squaredMagnitude(avg_tgt) < 1.0e6 * DBL_EPSILON * DBL_EPSILON) {
          avg_tgt = Normal(0.25 * tgt_0 + 0.75 * tgt_1);
        }
        avg_tgt.normalize();
        Pt proj_test =
            projectPtAlongHalfLineOntoParaboloid(a_paraboloid, avg_tgt, avg_pt);
        // Check if test projected point is inside the face
        UnsignedIndex_t best_proj_dir = 0;
        if (std::fabs(face_normal[best_proj_dir]) < std::fabs(face_normal[1]))
          best_proj_dir = 1;
        if (std::fabs(face_normal[best_proj_dir]) < std::fabs(face_normal[2]))
          best_proj_dir = 2;
        const UnsignedIndex_t start_id_for_search = (best_proj_dir + 1) % 3;
        auto current_half_edge = a_start;
        bool pt_internal_to_polygon = false;
        do {
          const Pt& location_0 =
              current_half_edge->getPreviousVertex()->getLocation().getPt();
          const Pt& location_1 =
              current_half_edge->getVertex()->getLocation().getPt();
          if (isPtBeforeIntersectionWithEdgeWithComponent(
                  proj_test, location_0, location_1, start_id_for_search)) {
            pt_internal_to_polygon = !pt_internal_to_polygon;
          }
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != a_start);

        // If test point is not inside face, invert tangents
        if (!pt_internal_to_polygon) {
          tgt_0 = -tgt_0;
          tgt_1 = -tgt_1;
        }
      }
    }
    return computeType3ContributionWithSplit<ReturnType>(
        a_paraboloid, face_normal, pt_1, pt_1, pt_0, tgt_1, tgt_0, a_surface);
  }

  return ReturnType::fromScalarConstant(0.0);
}  // namespace IRL

template <class ReturnType, class HalfEdgeType, class SurfaceOutputType>
ReturnType orientAndApplyType3CorrectionWithGradients(
    const AlignedParaboloid& a_paraboloid, HalfEdgeType* a_start,
    HalfEdgeType* a_end, SurfaceOutputType* a_surface) {
  using gradient_type = typename ReturnType::gradient_type;
  using pt_type = PtWithGradient<gradient_type>;
  const pt_type& pt_withgrad_0 = a_start->getVertex()->getLocation();
  const pt_type& pt_withgrad_1 = a_end->getVertex()->getLocation();
  const pt_type edge_vector_withgrad = pt_withgrad_1 - pt_withgrad_0;
  const Pt& pt_0 = pt_withgrad_0.getPt();
  const Pt& pt_1 = pt_withgrad_1.getPt();
  const Normal edge_vector = Normal(pt_1 - pt_0);
  const auto& face_plane = a_end->getFace()->getPlane();
  const auto& face_normal = face_plane.normal();

  // Compute plane gradient
  auto face_normal_withgrad =
      pt_type(Pt(face_normal[0], face_normal[1], face_normal[2]));
  auto& face_normal_grad = face_normal_withgrad.getData();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    face_normal_grad[d].setGradRx(0.0);
    face_normal_grad[d].setGradRy(0.0);
    face_normal_grad[d].setGradRz(0.0);
  }

  // Compute tangents and their gradient
  auto tgt_withgrad_0 = computeTangentVectorAndGradientAtPoint<pt_type>(
      a_paraboloid, face_normal_withgrad, pt_withgrad_0);
  auto tgt_withgrad_1 = computeTangentVectorAndGradientAtPoint<pt_type>(
      a_paraboloid, face_normal_withgrad, pt_withgrad_1);
  auto tgt_0 = Normal::fromPt(tgt_withgrad_0.getPt());
  auto tgt_1 = Normal::fromPt(tgt_withgrad_1.getPt());
  const bool elliptic_face =
      a_paraboloid.a() * a_paraboloid.b() > 0.0 && face_normal[2] != 0.0;

  if (!elliptic_face)  // The arc is from a hyperbola or parabola
  {
    const Normal n_cross_t0 = crossProduct(face_normal, tgt_0);
    if (std::fabs(n_cross_t0 * tgt_1) <
        100.0 * DBL_EPSILON)  // The tangents are (close to being) parallel
    {
      // We consider an arc with zero contribution
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        a_surface->addArc(
            RationalBezierArc(pt_1, 0.5 * (pt_0 + pt_1), pt_0, 0.0));
      }
      return ReturnType::fromScalarConstant(0.0);
    } else {  // The tangents are NOT parallel
      // Compute control point
      const double lambda_1 =
          -(n_cross_t0 * edge_vector) / safelyEpsilon(n_cross_t0 * tgt_1);
      const auto control_pt = Pt(pt_1 + lambda_1 * tgt_1);
      // Orient tangents to point towards control point
      tgt_withgrad_0 = (tgt_0 * Normal(control_pt - pt_0)) < 0.0
                           ? -tgt_withgrad_0
                           : tgt_withgrad_0;
      tgt_withgrad_1 = (tgt_1 * Normal(control_pt - pt_1)) < 0.0
                           ? -tgt_withgrad_1
                           : tgt_withgrad_1;
      // Construct Bezier arc from tangents
      const auto arc_with_gradient = RationalBezierArcWithGradient<pt_type>(
          pt_withgrad_1, tgt_withgrad_1, pt_withgrad_0, tgt_withgrad_0,
          face_normal_withgrad, a_paraboloid);
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        const auto arc =
            RationalBezierArc(pt_1, arc_with_gradient.control_point().getPt(),
                              pt_0, arc_with_gradient.weight());
        a_surface->addArc(arc);
      }
      return computeType3Contribution<ReturnType>(a_paraboloid,
                                                  arc_with_gradient);
    }
  } else  // The arc is from an ellipse
  {
    // Compute edge vectors and check if tangent is parallel to edge
    Normal edge_00 =
        Normal::fromPt(a_start->getVertex()->getLocation().getPt() -
                       a_start->getPreviousVertex()->getLocation().getPt());
    Normal edge_01 = Normal::fromPt(
        a_start->getNextHalfEdge()->getVertex()->getLocation().getPt() -
        a_start->getVertex()->getLocation().getPt());
    Normal edge_0 = (squaredMagnitude(edge_00) > squaredMagnitude(edge_01))
                        ? edge_00
                        : edge_01;
    edge_0.normalize();
    Normal edge_10 =
        Normal::fromPt(a_end->getVertex()->getLocation().getPt() -
                       a_end->getPreviousVertex()->getLocation().getPt());
    Normal edge_11 = Normal::fromPt(
        a_end->getNextHalfEdge()->getVertex()->getLocation().getPt() -
        a_end->getVertex()->getLocation().getPt());
    Normal edge_1 = (squaredMagnitude(edge_10) > squaredMagnitude(edge_11))
                        ? edge_10
                        : edge_11;
    edge_1.normalize();
    bool tgt_0_parallel_edge_0 =
        std::abs(1.0 - std::fabs(tgt_0 * edge_0)) < 100.0 * DBL_EPSILON;
    bool tgt_1_parallel_edge_1 =
        std::abs(1.0 - std::fabs(tgt_1 * edge_1)) < 100.0 * DBL_EPSILON;
    if (!tgt_0_parallel_edge_0 &&
        !tgt_1_parallel_edge_1)  // We orient the tangent with the edge
                                 // normal
    {
      const Normal edge_normal_0 = crossProduct(face_normal, edge_0);
      const Normal edge_normal_1 = crossProduct(face_normal, edge_1);
      tgt_withgrad_0 =
          (edge_normal_0 * tgt_0 < 0.0) ? -tgt_withgrad_0 : tgt_withgrad_0;
      tgt_withgrad_1 =
          (edge_normal_1 * tgt_1 < 0.0) ? -tgt_withgrad_1 : tgt_withgrad_1;
    } else {
      // Compute ellipse center an orient tangents accordingly
      const Pt conic_center = conicCenter(face_plane, a_paraboloid);
      auto center_to_pt_0 = Normal(pt_0 - conic_center);
      center_to_pt_0.normalize();
      auto center_to_pt_1 = Normal(pt_1 - conic_center);
      center_to_pt_1.normalize();
      Normal dummy_tgt_0 = crossProduct(face_normal, center_to_pt_0);
      Normal dummy_tgt_1 = crossProduct(face_normal, center_to_pt_1);
      assert(std::abs(tgt_0 * dummy_tgt_0) > 10.0 * DBL_EPSILON);
      assert(std::abs(tgt_1 * dummy_tgt_1) > 10.0 * DBL_EPSILON);
      tgt_0 = (tgt_0 * dummy_tgt_0) < 0.0 ? -tgt_0 : tgt_0;
      tgt_withgrad_0 =
          (tgt_0 * dummy_tgt_0) < 0.0 ? -tgt_withgrad_0 : tgt_withgrad_0;
      tgt_1 = (tgt_1 * dummy_tgt_1) > 0.0 ? -tgt_1 : tgt_1;
      tgt_withgrad_1 =
          (tgt_1 * dummy_tgt_1) > 0.0 ? -tgt_withgrad_1 : tgt_withgrad_1;
      // At this point, the tangents form a valid arc (but they may be
      // oriented in the wrong direction)
      if (!tgt_0_parallel_edge_0) {
        const Normal edge_normal_0 = crossProduct(face_normal, edge_0);
        if (edge_normal_0 * tgt_0 < 0.0) {
          tgt_withgrad_0 = -tgt_withgrad_0;
          tgt_withgrad_1 = -tgt_withgrad_1;
        }
      } else if (!tgt_1_parallel_edge_1) {
        const Normal edge_normal_1 = crossProduct(face_normal, edge_1);
        if (edge_normal_1 * tgt_1 < 0.0) {
          tgt_withgrad_0 = -tgt_withgrad_0;
          tgt_withgrad_1 = -tgt_withgrad_1;
        }
      } else {
        // The arcs should have been sorted in this case!! This allows to
        // project the mid-point onto the arc, and check that it belongs to
        // the face
        const Pt avg_pt = 0.5 * (pt_0 + pt_1);
        auto avg_tgt = Normal(0.5 * (tgt_0 + tgt_1));
        if (squaredMagnitude(avg_tgt) < 1.0e6 * DBL_EPSILON * DBL_EPSILON) {
          avg_tgt = Normal(0.25 * tgt_0 + 0.75 * tgt_1);
        }
        avg_tgt.normalize();
        Pt proj_test =
            projectPtAlongHalfLineOntoParaboloid(a_paraboloid, avg_tgt, avg_pt);
        // Check if test projected point is inside the face
        UnsignedIndex_t best_proj_dir = 0;
        if (std::fabs(face_normal[best_proj_dir]) < std::fabs(face_normal[1]))
          best_proj_dir = 1;
        if (std::fabs(face_normal[best_proj_dir]) < std::fabs(face_normal[2]))
          best_proj_dir = 2;
        const UnsignedIndex_t start_id_for_search = (best_proj_dir + 1) % 3;
        auto current_half_edge = a_start;
        bool pt_internal_to_polygon = false;
        do {
          const Pt& location_0 =
              current_half_edge->getPreviousVertex()->getLocation().getPt();
          const Pt& location_1 =
              current_half_edge->getVertex()->getLocation().getPt();
          if (isPtBeforeIntersectionWithEdgeWithComponent(
                  proj_test, location_0, location_1, start_id_for_search)) {
            pt_internal_to_polygon = !pt_internal_to_polygon;
          }
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != a_start);

        // If test point is not inside face, invert tangents
        if (!pt_internal_to_polygon) {
          tgt_withgrad_0 = -tgt_withgrad_0;
          tgt_withgrad_1 = -tgt_withgrad_1;
        }
      }
    }
    return computeType3ContributionWithGradientWithSplit<ReturnType>(
        a_paraboloid, face_normal_withgrad, pt_withgrad_1, pt_withgrad_1,
        pt_withgrad_0, tgt_withgrad_1, tgt_withgrad_0, a_surface);
  }

  return ReturnType::fromScalarConstant(0.0);
}  // namespace IRL

template <class HalfEdgeType, class PtType>
Normal determineNudgeDirection(const AlignedParaboloid& a_paraboloid,
                               HalfEdgeType* a_current_edge,
                               const PtType& a_inter_pt) {
  const auto& inter_pt = a_inter_pt.getPt();
  Normal nudge_direction = Normal(2.0 * a_paraboloid.a() * inter_pt[0],
                                  2.0 * a_paraboloid.b() * inter_pt[1], 1.0);
  Normal avg_face_normal =
      a_current_edge->getFace()->getPlane().normal() +
      a_current_edge->getOppositeHalfEdge()->getFace()->getPlane().normal();
  if (nudge_direction * avg_face_normal > 0.0) nudge_direction *= -1.0;
  nudge_direction.normalize();
  return nudge_direction;
}

// Note: This essentially abandons the half edges and vertices
// in a_complete_polytope, doing nothing to reclaim their memory usage.
// This should not be a problem unless resetPolyhedron is called
// many times for the same polyhedron intersection operation
template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, void>
resetPolyhedron(SegmentedHalfEdgePolyhedronType* a_polytope,
                HalfEdgePolytopeType* a_complete_polytope) {
  UnsignedIndex_t original_verts = 0;
  for (UnsignedIndex_t v = 0; v < a_polytope->getNumberOfVertices(); ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    if (vertex.needsToSeek()) {
      ++original_verts;
    } else {
      // New vertex from intersection, remove it and patch half-edges
      auto current_edge = vertex.getHalfEdge();
      auto original_edge = current_edge->getNextHalfEdge();
      doubleLinkHalfEdges(current_edge->getPreviousHalfEdge(), original_edge);
      original_edge->getFace()->setStartingHalfEdge(original_edge);
      auto opposite_edge = current_edge->getOppositeHalfEdge();
      doubleLinkHalfEdges(
          original_edge->getOppositeHalfEdge()->getPreviousHalfEdge(),
          opposite_edge);
      opposite_edge->getFace()->setStartingHalfEdge(opposite_edge);
      setMutualOpposites(original_edge, opposite_edge);
    }
  }
  a_polytope->setNumberOfVertices(original_verts);
  for (UnsignedIndex_t f = 0; f < a_polytope->getNumberOfFaces(); ++f) {
    (*a_polytope)[f]->clearIntersections();
  }
}

template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
void nudgePolyhedron(SegmentedHalfEdgePolyhedronType* a_polytope,
                     HalfEdgePolytopeType* a_complete_polytope,
                     const Normal a_nudge_direction,
                     const double a_nudge_epsilon,
                     const UnsignedIndex_t a_nudge_iter) {
  std::cout << " -----> Nudge along the direction " << a_nudge_direction
            << std::endl;
  using pt_type =
      typename SegmentedHalfEdgePolyhedronType::vertex_type::pt_type;
  for (UnsignedIndex_t n = 0; n < a_polytope->getNumberOfVertices(); ++n) {
    Pt shifted_location = a_polytope->getVertex(n)->getLocation().getPt() +
                          std::pow(1.1, static_cast<double>(a_nudge_iter)) *
                              a_nudge_epsilon * a_nudge_direction;
    a_polytope->getVertex(n)->setLocation(pt_type(shifted_location));
  }
}

// Assumes paraboloid of function 0 = a*x^2 + b*y^2 + z.
// We will truncate the polyhedron to exist in the region
// q < a*x^2 + b*y^2 + z
template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class SurfaceOutputType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
formParaboloidIntersectionBases(SegmentedHalfEdgePolyhedronType* a_polytope,
                                HalfEdgePolytopeType* a_complete_polytope,
                                const AlignedParaboloid& a_aligned_paraboloid,
                                const UnsignedIndex_t a_nudge_iter,
                                SurfaceOutputType* a_surface) {
  using vertex_type = typename SegmentedHalfEdgePolyhedronType::vertex_type;
  using pt_type = typename vertex_type::pt_type;
  using face_type = typename SegmentedHalfEdgePolyhedronType::face_type;
  using half_edge_type = typename HalfEdgePolytopeType::half_edge_type;

  assert(!(a_surface != nullptr &&
           std::is_same<SurfaceOutputType, NoSurfaceOutput>::value));

  ReturnType full_moments = ReturnType::fromScalarConstant(0.0);

  // Initialising variables for handling degenerate cases
  bool requires_nudge = false;
  const double nudge_epsilon = 10.0 * DBL_EPSILON;
  Normal nudge_direction;

  // Identify elliptic case
  const bool elliptic =
      a_aligned_paraboloid.a() * a_aligned_paraboloid.b() > 0.0;

  // Mark vertices as clipped/unclipped
  const auto starting_number_of_vertices = a_polytope->getNumberOfVertices();
  const auto starting_number_of_faces = a_polytope->getNumberOfFaces();
  UnsignedIndex_t number_of_vertices_above = 0;
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    vertex.setAsUnnecessaryToSeek();  // Reset all
    if (vertexBelow(vertex.getLocation().getPt(), a_aligned_paraboloid)) {
      vertex.markToBeNotClipped();
    } else {
      vertex.markToBeClipped();
      ++number_of_vertices_above;
    }
  }

  // Early termination cases, only possible with elliptic
  if (elliptic && a_aligned_paraboloid.a() > 0.0 &&
      number_of_vertices_above == 0) {
    // Whole volume below
    return ReturnType::calculateMoments(a_polytope);
  }

  if (elliptic && a_aligned_paraboloid.a() < 0.0 &&
      number_of_vertices_above == starting_number_of_vertices) {
    // Zero volume - will be current value of full_moments
    return full_moments;
  }

  // Clear visitation knowledge from polytope.
  for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
    (*a_polytope)[f]->markAsNotVisited();
    (*a_polytope)[f]->clearIntersections();
  }

  // Will have 0, 1, or 2 intercepts.
  // If intercept exists, place into HalfEdgeStructure
  const bool check_from_clipped = true;
  //      !elliptic || (elliptic && a_aligned_paraboloid.a() > 0.0);
  const bool check_from_unclipped = true;
  //      !elliptic || (elliptic && !check_from_clipped);

  auto edge_parallel_to_paraboloid =
      [&a_aligned_paraboloid, elliptic](
          const pt_type& a_edge_0, const pt_type& a_edge_1,
          const pt_type& a_intersection) -> UnsignedIndex_t {
    if (!elliptic) {
      return 0;
    }
    const Pt line = a_edge_1.getPt() - a_edge_0.getPt();
    const auto edge_direction = Normal::fromPtNormalized(line);
    const auto paraboloid_normal = getParaboloidSurfaceNormal(
        a_aligned_paraboloid, a_intersection.getPt());
    return std::fabs(edge_direction * paraboloid_normal) < 1000.0 * DBL_EPSILON
               ? 1
               : 0;
  };

  auto edge_parallel_to_paraboloid2 =
      [&a_aligned_paraboloid, elliptic](
          const pt_type& a_edge_0, const pt_type& a_edge_1,
          const pt_type& a_intersection_0,
          const pt_type& a_intersection_1) -> UnsignedIndex_t {
    if (!elliptic) {
      return 0;
    }
    UnsignedIndex_t parallel = 0;

    const Pt line = a_edge_1.getPt() - a_edge_0.getPt();
    const auto edge_direction = Normal::fromPtNormalized(line);
    auto paraboloid_normal = getParaboloidSurfaceNormal(
        a_aligned_paraboloid, a_intersection_0.getPt());
    if (std::fabs(edge_direction * paraboloid_normal) < 1000.0 * DBL_EPSILON) {
      ++parallel;
    }
    paraboloid_normal = getParaboloidSurfaceNormal(a_aligned_paraboloid,
                                                   a_intersection_1.getPt());
    if (std::fabs(edge_direction * paraboloid_normal) < 1000.0 * DBL_EPSILON) {
      ++parallel;
    }
    return parallel;
  };

  // Compute gradients of polyhedron corners (if gradients are requested)
  if constexpr (has_embedded_gradient<ReturnType>::value) {
    for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
      auto& gradient = a_polytope->getVertex(v)->getLocation().getData();
      gradient[0].setGradTx(-1.0);
      gradient[1].setGradTy(-1.0);
      gradient[2].setGradTz(-1.0);
    }
  }

  // Compute intersections
  StackVector<std::pair<pt_type, double>, 2> edge_intercepts;
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    vertex.setToSeek();
    if (check_from_clipped && vertex.isClipped()) {
      // CASE WHERE STARTING VERTEX IS CLIPPED
      auto current_edge = vertex.getHalfEdge();
      const auto starting_edge = current_edge;
      do {
        // If it has needsToSeek set, it means it is a new vertex or
        // already visited. Either way, do not need to check for
        // intersection
        if (current_edge->getPreviousVertex()->needsToSeek()) {
          current_edge =
              current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
          continue;
        }
        // If previous vertex is not clipped, single-intercept
        // If previous vertex is clipped, need to check for double-intercept
        // Checking for double-intercept and calculating single-intercept
        // is the same routine, so just always do it.
        const auto& edge_start =
            current_edge->getPreviousVertex()->getLocation();
        const auto& edge_end = current_edge->getVertex()->getLocation();
        checkAndFindIntercepts<pt_type>(a_aligned_paraboloid, edge_start,
                                        edge_end, &edge_intercepts,
                                        nudge_epsilon);

        // assert(current_edge->getPreviousVertex()->isClipped() ||
        //        edge_intercepts.size() == 1);
        // Size of returned intercepts indicates single or double intercept
        // (or none)
        if (edge_intercepts.size() == 1) {
          // Check for intersection near end point
          if (edge_intercepts[0].second < nudge_epsilon ||
              (1.0 - edge_intercepts[0].second) < nudge_epsilon) {
            requires_nudge = true;
            nudge_direction = determineNudgeDirection(
                a_aligned_paraboloid, current_edge, edge_intercepts[0].first);
            break;
          }

          // Add vertex to half edge structure, resetting connectivity and
          // creating a new half edge (and new opposite half edge)
          placeSingleIntercept(edge_intercepts[0].first, current_edge,
                               a_polytope, a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->addIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->setStartingHalfEdge(opposite_half_edge);
          opposite_face->addIntersection();

          auto edge_parallel_intersections = edge_parallel_to_paraboloid(
              edge_start, edge_end, edge_intercepts[0].first);
          if (edge_parallel_intersections > 0) {
            current_face->addEdgeParallelIntersections(
                edge_parallel_intersections);
            opposite_face->addEdgeParallelIntersections(
                edge_parallel_intersections);
          }

        } else if (edge_intercepts.size() == 2) {
          // Check for intersection near end point
          if (edge_intercepts[0].second < nudge_epsilon) {
            requires_nudge = true;
            nudge_direction = determineNudgeDirection(
                a_aligned_paraboloid, current_edge, edge_intercepts[0].first);
            break;
          }
          if ((1.0 - edge_intercepts[1].second) < nudge_epsilon) {
            requires_nudge = true;
            nudge_direction = determineNudgeDirection(
                a_aligned_paraboloid, current_edge, edge_intercepts[1].first);
            break;
          }

          placeDoubleIntercept(edge_intercepts, current_edge, a_polytope,
                               a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->setStartingHalfEdge(
              current_edge->getPreviousHalfEdge()->getPreviousHalfEdge());
          current_face->addDoubleIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->setStartingHalfEdge(
              current_edge->getOppositeHalfEdge());
          opposite_face->addDoubleIntersection();

          auto edge_parallel_intersections = edge_parallel_to_paraboloid2(
              edge_start, edge_end, edge_intercepts[0].first,
              edge_intercepts[1].first);
          if (edge_parallel_intersections > 0) {
            current_face->addEdgeParallelIntersections(
                edge_parallel_intersections);
            opposite_face->addEdgeParallelIntersections(
                edge_parallel_intersections);
          }
        }
        current_edge =
            current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
      } while (current_edge != starting_edge);
    } else if (check_from_unclipped && vertex.isNotClipped()) {
      // CASE WHERE STARTING VERTEX IS UNCLIPPED
      auto current_edge = vertex.getHalfEdge();
      const auto starting_edge = current_edge;
      do {
        // If it has needsToSeek set, it means it is a new vertex or
        // already visited. Either way, do not need to check for
        // intersection
        if (current_edge->getPreviousVertex()->needsToSeek()) {
          current_edge =
              current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
          continue;
        }
        // If previous vertex is clipped, single-intercept
        // If previous vertex is not clipped, need to check for
        // double-intercept Checking for double-intercept and calculating
        // single-intercept is the same routine, so just always do it.
        const auto& edge_start =
            current_edge->getPreviousVertex()->getLocation();
        const auto& edge_end = current_edge->getVertex()->getLocation();
        checkAndFindIntercepts<pt_type>(a_aligned_paraboloid, edge_start,
                                        edge_end, &edge_intercepts,
                                        nudge_epsilon);

        // assert(current_edge->getPreviousVertex()->isNotClipped() ||
        //        edge_intercepts.size() == 1);
        // Size of returned intercepts indicates single or double intercept
        // (or none)
        if (edge_intercepts.size() == 1) {
          // Check for intersection near end point
          if (edge_intercepts[0].second < nudge_epsilon ||
              (1.0 - edge_intercepts[0].second) < nudge_epsilon) {
            requires_nudge = true;
            nudge_direction = determineNudgeDirection(
                a_aligned_paraboloid, current_edge, edge_intercepts[0].first);
            break;
          }

          // Add vertex to half edge structure, resetting connectivity and
          // creating a new half edge (and new opposite half edge)
          placeSingleIntercept(edge_intercepts[0].first, current_edge,
                               a_polytope, a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->setStartingHalfEdge(
              current_edge->getPreviousHalfEdge());
          current_face->addIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->addIntersection();

          auto edge_parallel_intersections = edge_parallel_to_paraboloid(
              edge_start, edge_end, edge_intercepts[0].first);
          if (edge_parallel_intersections > 0) {
            current_face->addEdgeParallelIntersections(
                edge_parallel_intersections);
            opposite_face->addEdgeParallelIntersections(
                edge_parallel_intersections);
          }

        } else if (edge_intercepts.size() == 2) {
          // Check for intersection near end point
          if (edge_intercepts[0].second < nudge_epsilon) {
            requires_nudge = true;
            nudge_direction = determineNudgeDirection(
                a_aligned_paraboloid, current_edge, edge_intercepts[0].first);
            break;
          }
          if ((1.0 - edge_intercepts[1].second) < nudge_epsilon) {
            requires_nudge = true;
            nudge_direction = determineNudgeDirection(
                a_aligned_paraboloid, current_edge, edge_intercepts[1].first);
            break;
          }

          placeDoubleIntercept(edge_intercepts, current_edge, a_polytope,
                               a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->setStartingHalfEdge(
              current_edge->getPreviousHalfEdge());
          current_face->addDoubleIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->setStartingHalfEdge(
              opposite_half_edge->getNextHalfEdge());
          opposite_face->addDoubleIntersection();

          auto edge_parallel_intersections = edge_parallel_to_paraboloid2(
              edge_start, edge_end, edge_intercepts[0].first,
              edge_intercepts[1].first);
          if (edge_parallel_intersections > 0) {
            current_face->addEdgeParallelIntersections(
                edge_parallel_intersections);
            opposite_face->addEdgeParallelIntersections(
                edge_parallel_intersections);
          }
        }
        current_edge =
            current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
      } while (current_edge != starting_edge);
    }
    if (requires_nudge) {
      for (UnsignedIndex_t i = v + 1; i < starting_number_of_vertices; ++i) {
        a_polytope->getVertex(i)->setToSeek();
      }
      break;
    }
  }
  assert(a_polytope->checkValidHalfEdgeStructure());
  // After leaving above loop, all faces that were intersected have a
  // starting half edge that ends on a new (intersection) vertex and the
  // previous vertex is unclipped (or new). (i.e., the edge exists in the
  // unclipped portion.
  // All new vertices also have their half edge being one that has its
  // edge in the unclipped portion.

  const auto vertices_after_intersection = a_polytope->getNumberOfVertices();
  const auto new_intersection_vertices =
      vertices_after_intersection - starting_number_of_vertices;

  for (UnsignedIndex_t v = starting_number_of_vertices;
       v < vertices_after_intersection; ++v) {
    // Original vertices will be set as needsToSeek()
    // Reset newly created vertices will have doesNotNeedToSeek()
    a_polytope->getVertex(v)->setAsUnnecessaryToSeek();
  }

  // If intersection is too close to corners, shift polyhedron and try again
  if (requires_nudge) {
    // Clean half-edge structure by removing intersections
    if (new_intersection_vertices > 0) {
      resetPolyhedron(a_polytope, a_complete_polytope);
    }
    nudgePolyhedron(a_polytope, a_complete_polytope, nudge_direction,
                    nudge_epsilon, a_nudge_iter + 1);
    return formParaboloidIntersectionBases<ReturnType>(
        a_polytope, a_complete_polytope, a_aligned_paraboloid, a_nudge_iter + 1,
        a_surface);
  }

  // Check for face-only intersections. Can only happen with elliptic
  if (elliptic) {
    if (a_aligned_paraboloid.a() > 0.0) {
      for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
        auto& face = *(*a_polytope)[f];
        if (face.getNumberOfIntersections() > 0) {
          continue;
        }
        // Face will be FaceOnly intersect if
        // 1 - Any vertex on the face is below the paraboloid
        // 2 - No z-component on the face is below the maximum for
        // a downward opening elliptic paraboloid, which is 0.
        bool face_valid = false;
        const auto starting_half_edge = face.getStartingHalfEdge();
        auto current_half_edge = starting_half_edge;
        do {
          const auto& vertex = *(current_half_edge->getVertex());
          if (vertex.isNotClipped()) {
            face_valid = false;
            break;
          }
          if (vertex.getLocation().getPt()[2] < 0.0) {
            face_valid = true;
          }
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != starting_half_edge);
        if (!face_valid) {
          continue;
        }

        // Made it this far, check if there is a face-only intersection
        const auto& face_plane = face.getPlane();
        if (std::fabs(face_plane.normal()[2]) > DBL_EPSILON) {
          // Get ellipse on this face
          if (ellipseContainedInFace(a_aligned_paraboloid, face_plane,
                                     starting_half_edge)) {
            full_moments += computeFaceOnlyContribution<ReturnType>(
                a_aligned_paraboloid, face_plane,
                starting_half_edge->getVertex()->getLocation());
            // Return surface parametrization
            if constexpr (!std::is_same<SurfaceOutputType,
                                        NoSurfaceOutput>::value) {
              addEllipseToSurfaceOutput(a_aligned_paraboloid, face_plane,
                                        a_surface);
            }
          }
        }
      }
    } else {
      // Case where a_aligned_paraboloid.a() < 0.0)
      for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
        auto& face = *(*a_polytope)[f];
        if (face.getNumberOfIntersections() > 0) {
          continue;
        }
        // Face will be FaceOnly intersect if
        // 1 - Any vertex on the face is above the paraboloid
        // 2 - No z-component on the face is above the minimum for
        // an upward opening elliptic paraboloid, which is 0.
        bool face_valid = false;
        const auto starting_half_edge = face.getStartingHalfEdge();
        auto current_half_edge = starting_half_edge;
        do {
          const auto& vertex = *(current_half_edge->getVertex());
          if (vertex.isClipped()) {
            face_valid = false;
            break;
          }
          if (vertex.getLocation().getPt()[2] > 0.0) {
            face_valid = true;
          }
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != starting_half_edge);
        if (!face_valid) {
          continue;
        }

        // Made it this far, check if there is a face-only intersection
        const auto& face_plane = face.getPlane();
        if (std::fabs(face_plane.normal()[2]) > DBL_EPSILON) {
          // Get ellipse on this face
          if (ellipseContainedInFace(a_aligned_paraboloid, face_plane,
                                     starting_half_edge)) {
            full_moments += computeFaceOnlyContribution<ReturnType>(
                a_aligned_paraboloid, face_plane,
                starting_half_edge->getVertex()->getLocation());
            if constexpr (!std::is_same<SurfaceOutputType,
                                        NoSurfaceOutput>::value) {
              addEllipseToSurfaceOutput(a_aligned_paraboloid, face_plane,
                                        a_surface);
            }
          }
        }
      }
    }
  }

  if (new_intersection_vertices == 0) {
    if (number_of_vertices_above == starting_number_of_vertices) {
      // All points above
      return full_moments;
    } else {
      // All points below
      return ReturnType::calculateMoments(a_polytope) + full_moments;
    }
  }

  for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
    auto& face = *(*a_polytope)[f];
    const auto& face_normal = face.getPlane().normal();
    auto starting_half_edge = face.getStartingHalfEdge();

    const auto intersection_size = face.getNumberOfIntersections();
    if (intersection_size == 0) {
      if (starting_half_edge->getVertex()->isNotClipped()) {
        // The face is entirely below
        const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
        auto current_half_edge =
            starting_half_edge->getNextHalfEdge()->getNextHalfEdge();
        auto prev_pt = current_half_edge->getPreviousVertex()->getLocation();
        do {
          const auto& curr_pt = current_half_edge->getVertex()->getLocation();
          full_moments +=
              computeType1Contribution<ReturnType>(ref_pt, prev_pt, curr_pt);
          prev_pt = curr_pt;
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != starting_half_edge);
      }
    } else if (intersection_size == 2) {
      const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
      half_edge_type* exit_half_edge;
      full_moments += computeUnclippedSegmentType1Contribution<ReturnType>(
          a_aligned_paraboloid, ref_pt, starting_half_edge, exit_half_edge,
          true);
      full_moments += computeNewEdgeSegmentContribution<ReturnType>(
          a_aligned_paraboloid, ref_pt, starting_half_edge, exit_half_edge,
          true, a_surface);
    } else {
      // More than 2 intersections. More complicated cases.
      const bool elliptic_face =
          a_aligned_paraboloid.a() * a_aligned_paraboloid.b() > 0.0 &&
          face_normal[2] != 0.0;
      const bool hyperbolic_face =
          a_aligned_paraboloid.a() * a_aligned_paraboloid.b() < 0.0 &&
          face_normal[2] != 0.0;

      // -----> If NO surface ouput : DOES NOT NEED SORTING
      if (face.getNumberOfEdgeParallelIntersections() < 2 &&
          std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        const bool reverse = a_aligned_paraboloid.a() < 0.0;
        if (elliptic_face) {
          const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
          half_edge_type* exit_half_edge;
          auto current_edge = starting_half_edge;
          bool skip_first = true;
          std::size_t found_intersections = 0;
          do {
            full_moments +=
                computeUnclippedSegmentType1Contribution<ReturnType>(
                    a_aligned_paraboloid, ref_pt, current_edge, exit_half_edge,
                    skip_first);

            if (reverse) {
              full_moments += computeNewEdgeSegmentContribution<ReturnType>(
                  a_aligned_paraboloid, ref_pt, current_edge, exit_half_edge,
                  skip_first, a_surface);
              current_edge = exit_half_edge->getNextHalfEdge();
              while (current_edge->getVertex()->needsToSeek()) {
                current_edge = current_edge->getNextHalfEdge();
              }
            } else {
              current_edge = exit_half_edge->getNextHalfEdge();
              while (current_edge->getVertex()->needsToSeek()) {
                current_edge = current_edge->getNextHalfEdge();
              }
              full_moments += computeNewEdgeSegmentContribution<ReturnType>(
                  a_aligned_paraboloid, ref_pt, current_edge, exit_half_edge,
                  current_edge == starting_half_edge, a_surface);
            }
            skip_first = false;
            found_intersections += 2;
          } while (found_intersections != intersection_size);
        } else {
          SmallVector<half_edge_type*, 10> intersections;
          // Find intersections and determine status
          auto current_edge = starting_half_edge;
          bool reverse_order = false;
          half_edge_type* exit_half_edge;
          std::size_t found_intersections = 0;
          const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
          bool skip_first = true;
          do {
            full_moments +=
                computeUnclippedSegmentType1Contribution<ReturnType>(
                    a_aligned_paraboloid, ref_pt, current_edge, exit_half_edge,
                    skip_first);
            current_edge->getVertex()->markAsEntry();
            exit_half_edge->getVertex()->markAsExit();
            intersections.push_back(current_edge);
            intersections.push_back(exit_half_edge);
            current_edge = exit_half_edge->getNextHalfEdge();
            while (current_edge->getVertex()->needsToSeek()) {
              current_edge = current_edge->getNextHalfEdge();
            }
            found_intersections += 2;
            skip_first = false;
          } while (found_intersections != intersection_size);

          if (hyperbolic_face) {
            const std::array<double, 2> conic_center{
                {face_normal[0] /
                     (2.0 * a_aligned_paraboloid.a() * face_normal[2]),
                 face_normal[1] /
                     (2.0 * a_aligned_paraboloid.b() * face_normal[2])}};
            const auto& pt_in_plane =
                intersections[0]->getVertex()->getLocation().getPt();
            const double delta_face = (face_normal[0] * pt_in_plane[0] +
                                       face_normal[1] * pt_in_plane[1] +
                                       face_normal[2] * pt_in_plane[2]) /
                                      face_normal[2];
            const double gamma_face =
                a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
                a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] -
                delta_face;
            const double z_diff =
                face_normal[0] * conic_center[0] / face_normal[2] +
                face_normal[1] * conic_center[1] / face_normal[2] - gamma_face -
                2.0 * delta_face;
            const std::size_t split_ind =
                a_aligned_paraboloid.a() * gamma_face > 0.0 ? 0 : 1;

            bool vertices_on_same_branch = true;
            for (std::size_t i = 1; i < intersection_size; ++i) {
              const auto& curr_pt =
                  intersections[i]->getVertex()->getLocation().getPt();
              if ((pt_in_plane[split_ind] - conic_center[split_ind]) *
                      (curr_pt[split_ind] - conic_center[split_ind]) <
                  0.0) {
                vertices_on_same_branch = false;
                break;
              }
            }
            bool total_invert = z_diff < 0.0 ? !vertices_on_same_branch
                                             : vertices_on_same_branch;
            if (total_invert) {
              reverse_order = true;
            }
          } else {
            Pt intersection_avg = Pt(0.0, 0.0, 0.0);
            for (std::size_t i = 0; i < intersection_size; ++i) {
              const auto& curr_pt =
                  intersections[i]->getVertex()->getLocation().getPt();
              intersection_avg += curr_pt;
            }
            intersection_avg /= static_cast<double>(intersection_size);
            const double z_diff = intersection_avg[2] +
                                  a_aligned_paraboloid.a() *
                                      intersection_avg[0] *
                                      intersection_avg[0] +
                                  a_aligned_paraboloid.b() *
                                      intersection_avg[1] * intersection_avg[1];
            if (z_diff > 0.0) {
              reverse_order = true;
            }
          }

          // Traverse face from entry->exit
          // Integrate internal portion of face
          // Integrate wedge on new edge from exit->entry
          bool first_entry = true;
          // Loop over vertices, only use entry vertices
          for (std::size_t i = 0; i < intersection_size; i += 2) {
            // Identify entry vertex
            auto entry_half_edge = intersections[i];
            // Loop over internal portion from entry->exit
            int exit_index;
            if (reverse_order) {
              exit_index = i != intersection_size - 1 ? i + 1 : 0;
            } else {
              exit_index = i != 0 ? i - 1 : intersection_size - 1;
            }

            auto exit_half_edge = intersections[exit_index];
            full_moments += computeNewEdgeSegmentContribution<ReturnType>(
                a_aligned_paraboloid, ref_pt, entry_half_edge, exit_half_edge,
                first_entry, a_surface);
            first_entry = false;
          }
          // Clear list of intersections
          intersections.clear();
        }
      }
      // -----> If surface ouput : NEEDS SORTING
      else {
        using stype = std::pair<half_edge_type*, double>;
        SmallVector<stype, 10> intersections;
        // Find intersections and determine status
        auto current_edge = starting_half_edge;
        bool reverse_order = false;
        half_edge_type* exit_half_edge;
        std::size_t found_intersections = 0;
        const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
        bool skip_first = true;
        do {
          full_moments += computeUnclippedSegmentType1Contribution<ReturnType>(
              a_aligned_paraboloid, ref_pt, current_edge, exit_half_edge,
              skip_first);
          current_edge->getVertex()->markAsEntry();
          exit_half_edge->getVertex()->markAsExit();
          intersections.push_back(
              std::pair<half_edge_type*, double>({current_edge, DBL_MAX}));
          intersections.push_back(
              std::pair<half_edge_type*, double>({exit_half_edge, DBL_MAX}));
          current_edge = exit_half_edge->getNextHalfEdge();
          while (current_edge->getVertex()->needsToSeek()) {
            current_edge = current_edge->getNextHalfEdge();
          }
          found_intersections += 2;
          skip_first = false;
        } while (found_intersections != intersection_size);

        if (elliptic_face) {
          const std::array<double, 2> conic_center{
              {face_normal[0] /
                   (2.0 * a_aligned_paraboloid.a() * face_normal[2]),
               face_normal[1] /
                   (2.0 * a_aligned_paraboloid.b() * face_normal[2])}};
          const double normal_invert = std::copysign(1.0, face_normal[2]);
          const double invert =
              a_aligned_paraboloid.a() < 0.0 ? -normal_invert : normal_invert;

          for (auto& element : intersections) {
            const auto& pt = element.first->getVertex()->getLocation().getPt();
            element.second = invert * std::atan2(pt[1] - conic_center[1],
                                                 pt[0] - conic_center[0]);
          }
          std::sort(intersections.begin(), intersections.end(),
                    [](const stype& a, const stype& b) {
                      return a.second < b.second;
                    });
        } else if (hyperbolic_face) {
          std::size_t pos_end = 0;
          std::size_t neg_end = 0;
          auto intersection_copy = intersections;
          const std::array<double, 2> conic_center{
              {face_normal[0] /
                   (2.0 * a_aligned_paraboloid.a() * face_normal[2]),
               face_normal[1] /
                   (2.0 * a_aligned_paraboloid.b() * face_normal[2])}};
          const auto& pt_in_plane =
              intersections[0].first->getVertex()->getLocation().getPt();
          const double delta_face = (face_normal[0] * pt_in_plane[0] +
                                     face_normal[1] * pt_in_plane[1] +
                                     face_normal[2] * pt_in_plane[2]) /
                                    face_normal[2];
          const double gamma_face =
              a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
              a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] -
              delta_face;
          const std::size_t split_ind =
              a_aligned_paraboloid.a() * gamma_face > 0.0 ? 0 : 1;
          const std::size_t store_ind = split_ind == 0 ? 1 : 0;
          const double z_center_plane =
              -face_normal[0] * conic_center[0] / face_normal[2] -
              face_normal[1] * conic_center[1] / face_normal[2] + delta_face;
          const double z_center_paraboloid =
              -a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] -
              a_aligned_paraboloid.b() * conic_center[1] * conic_center[1];
          double total_invert = std::copysign(1.0, face_normal[2]);
          total_invert = z_center_plane < z_center_paraboloid ? total_invert
                                                              : -total_invert;
          total_invert = split_ind == 0 ? total_invert : -total_invert;
          for (auto& element : intersections) {
            const auto& pt = element.first->getVertex()->getLocation().getPt();
            if (pt[split_ind] > conic_center[split_ind]) {
              const auto loc = pos_end++;
              intersection_copy[loc].first = element.first;
              intersection_copy[loc].second = total_invert * pt[store_ind];
            } else {
              const auto loc = intersection_size - 1 - (neg_end++);
              intersection_copy[loc].first = element.first;
              intersection_copy[loc].second = -total_invert * pt[store_ind];
            }
          }
          assert(pos_end + neg_end == intersection_size);
          std::sort(intersection_copy.begin(),
                    intersection_copy.begin() + pos_end,
                    [](const stype& a, const stype& b) {
                      return a.second < b.second;
                    });
          std::sort(intersection_copy.begin() + pos_end,
                    intersection_copy.end(),
                    [](const stype& a, const stype& b) {
                      return a.second < b.second;
                    });
          intersections = intersection_copy;
        } else {
          Pt intersection_avg = Pt(0.0, 0.0, 0.0);
          for (std::size_t i = 0; i < intersection_size; ++i) {
            const auto& curr_pt =
                intersections[i].first->getVertex()->getLocation().getPt();
            intersection_avg += curr_pt;
          }
          intersection_avg /= static_cast<double>(intersection_size);
          const double z_diff = intersection_avg[2] +
                                a_aligned_paraboloid.a() * intersection_avg[0] *
                                    intersection_avg[0] +
                                a_aligned_paraboloid.b() * intersection_avg[1] *
                                    intersection_avg[1];
          if (face_normal[2] != 0.0) {
            const double normal_invert = std::copysign(1.0, face_normal[2]);
            const double invert = z_diff < 0.0 ? normal_invert : -normal_invert;
            for (auto& element : intersections) {
              const auto& pt =
                  element.first->getVertex()->getLocation().getPt();
              element.second = invert * std::atan2(pt[1] - intersection_avg[1],
                                                   pt[0] - intersection_avg[0]);
            }
            std::sort(intersections.begin(), intersections.end(),
                      [](const stype& a, const stype& b) {
                        return a.second < b.second;
                      });
          } else {
          }
        }
        // Traverse face from entry->exit
        // Integrate internal portion of face
        // Integrate wedge on new edge from exit->entry
        // Loop over vertices, only use entry vertices
        auto prev_vertex = intersections[0].first->getPreviousVertex();
        auto entry_half_edge = intersections[0].first;
        const bool entry_first =
            (prev_vertex->isClipped() ||
             (prev_vertex->doesNotNeedToSeek() &&
              entry_half_edge->getNextHalfEdge()->getVertex()->isNotClipped()));
        std::size_t start_id = entry_first ? 0 : 1;
        for (std::size_t i = start_id; i < intersection_size; i += 2) {
          // Identify entry vertex
          const auto entry_half_edge = intersections[i].first;
          // Loop over internal portion from entry->exit
          const UnsignedIndex_t exit_index =
              i != 0 ? i - 1 : intersection_size - 1;
          const auto exit_half_edge = intersections[exit_index].first;
          full_moments += computeNewEdgeSegmentContribution<ReturnType>(
              a_aligned_paraboloid, ref_pt, entry_half_edge, exit_half_edge,
              false, a_surface);
        }
        // Clear list of intersections
        intersections.clear();
      }
    }
  }

  return full_moments;
}
}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_TPP_
