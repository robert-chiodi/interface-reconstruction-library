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
#include <random>
#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_initializer.tpp"
#include "irl/generic_cutting/paraboloid_intersection/moment_contributions.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/geometry/half_edge_structures/brep_to_half_edge.h"
#include "irl/helpers/mymath.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/paraboloid_reconstruction/rational_bezier_arc.h"

namespace IRL {

template <class ScalarType>
inline const ScalarType distance_epsilon(void);

template <>
inline const double distance_epsilon(void) {
  return 1.0e2 * DBL_EPSILON;
}

template <>
inline const Quad_t distance_epsilon(void) {
  return 1.0e6q * FLT128_EPSILON;
}

template <class ScalarType>
inline const ScalarType angle_epsilon(void);

template <>
inline const double angle_epsilon(void) {
  return 1.0e6 * DBL_EPSILON;
}

template <>
inline const Quad_t angle_epsilon(void) {
  return 1.0e10q * FLT128_EPSILON;
}

/******************** Tangent to surface at given point ********************/
template <class ScalarType>
inline NormalBase<ScalarType> computeTangentVectorAtPoint(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const NormalBase<ScalarType>& a_plane_normal,
    const PtBase<ScalarType>& a_pt) {
  /* Defining constants */
  const ScalarType EPSILON = distance_epsilon<ScalarType>();
  const ScalarType ZERO = static_cast<ScalarType>(0);

  /* Function */
  NormalBase<ScalarType> surface_normal =
      getParaboloidSurfaceNormal(a_paraboloid, a_pt);
  surface_normal.approximatelyNormalize();
  NormalBase<ScalarType> tangent_at_pt =
      crossProduct(a_plane_normal, surface_normal);
  if (squaredMagnitude(tangent_at_pt) < EPSILON * EPSILON) {
    return NormalBase<ScalarType>(ZERO, ZERO, ZERO);
  }
  const ScalarType normal_correction = tangent_at_pt * a_plane_normal;
  tangent_at_pt = tangent_at_pt - normal_correction * a_plane_normal;
  // tangent_at_pt.normalize();
  return tangent_at_pt;
}

/************** Tangent to surface + gradient at given point **************/
template <class ScalarType, class PtWithGradientType>
inline PtWithGradientType computeTangentVectorAndGradientAtPoint(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtWithGradientType& a_plane_normal, const PtWithGradientType& a_pt) {
  /* Defining constants and types */
  using gradient_type = typename PtWithGradientType::gradient_type;
  const ScalarType EPSILON = distance_epsilon<ScalarType>();

  /* Function */
  const auto surface_normal_withgrad =
      getParaboloidSurfaceNormalWithGradient<PtWithGradientType>(a_paraboloid,
                                                                 a_pt);
  const NormalBase<ScalarType> surface_normal =
      NormalBase<ScalarType>::fromPt(surface_normal_withgrad.getPt());
  const auto& surface_normal_grad = surface_normal_withgrad.getData();
  const auto& a_plane_normal_grad = a_plane_normal.getData();
  const ScalarType tangent_at_pt_x = a_plane_normal[1] * surface_normal[2] -
                                     a_plane_normal[2] * surface_normal[1];
  const ScalarType tangent_at_pt_y = a_plane_normal[2] * surface_normal[0] -
                                     a_plane_normal[0] * surface_normal[2];
  const ScalarType tangent_at_pt_z = a_plane_normal[0] * surface_normal[1] -
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
  const ScalarType norm_tangent = sqrt(tangent_at_pt_x * tangent_at_pt_x +
                                       tangent_at_pt_y * tangent_at_pt_y +
                                       tangent_at_pt_z * tangent_at_pt_z);
  const auto norm_tangent_grad = (tangent_at_pt_gradx * tangent_at_pt_x +
                                  tangent_at_pt_grady * tangent_at_pt_y +
                                  tangent_at_pt_gradz * tangent_at_pt_z) /
                                 safelyEpsilon(norm_tangent);
  assert(norm_tangent > EPSILON);
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

/************** Tangent to surface + orient at given point **************/
template <class ScalarType>
inline NormalBase<ScalarType> computeAndCorrectTangentVectorAtPt(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const NormalBase<ScalarType>& a_plane_normal,
    const PtBase<ScalarType>& a_origin_pt, const PtBase<ScalarType>& a_end_pt,
    const NormalBase<ScalarType>& a_end_tangent,
    const PtBase<ScalarType>& a_intersection_pt) {
  /* Defining constants and types */
  const ScalarType ZERO = static_cast<ScalarType>(0);

  /* Function */
  NormalBase<ScalarType> tangent = computeTangentVectorAtPoint<ScalarType>(
      a_paraboloid, a_plane_normal, a_intersection_pt);
  tangent.normalize();
  const NormalBase<ScalarType> edge_normal =
      crossProduct(a_plane_normal, a_end_pt - a_intersection_pt);
  if ((a_end_tangent * edge_normal > ZERO) == (tangent * edge_normal > ZERO)) {
    tangent = -tangent;
  }
  return tangent;
}

/********** Tangent to surface + orient + gradient at given point **********/
template <class ScalarType, class PtTypeWithGradient>
inline PtTypeWithGradient computeAndCorrectTangentVectorAndGradientAtPt(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtTypeWithGradient& a_plane_normal,
    const PtTypeWithGradient& a_origin_pt, const PtTypeWithGradient& a_end_pt,
    const PtTypeWithGradient& a_end_tangent,
    const PtTypeWithGradient& a_intersection_pt) {
  /* Defining constants and types */
  const ScalarType ZERO = static_cast<ScalarType>(0);

  /* Function */
  PtTypeWithGradient tangent_withgrad =
      computeTangentVectorAndGradientAtPoint<ScalarType, PtTypeWithGradient>(
          a_paraboloid, a_plane_normal, a_intersection_pt);
  const NormalBase<ScalarType> tangent =
      NormalBase<ScalarType>::fromPt(tangent_withgrad.getPt());
  const NormalBase<ScalarType> plane_normal =
      NormalBase<ScalarType>::fromPt(a_plane_normal.getPt());
  const NormalBase<ScalarType> edge_vector1 =
      NormalBase<ScalarType>::fromPt(a_end_pt.getPt() - a_origin_pt.getPt());
  const NormalBase<ScalarType> edge_vector2 = NormalBase<ScalarType>::fromPt(
      a_intersection_pt.getPt() - a_origin_pt.getPt());
  const NormalBase<ScalarType> end_tangent =
      NormalBase<ScalarType>::fromPt(a_end_tangent.getPt());
  const NormalBase<ScalarType> correct_sign =
      crossProduct(plane_normal, edge_vector1);
  if ((end_tangent * correct_sign > ZERO) !=
      (tangent * Normal(crossProduct(plane_normal, edge_vector2)) > ZERO)) {
    tangent_withgrad = -tangent_withgrad;
  }
  return tangent_withgrad;
}

/***************** Calculate arc contribution (with split) *******************/
template <class ReturnType, class ScalarType, class SurfaceOutputType,
          class PtType>
ReturnType computeType3ContributionWithSplit(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const NormalBase<ScalarType>& a_plane_normal, const PtType& a_pt_ref,
    const PtType& a_pt_0, const PtType& a_pt_1,
    const NormalBase<ScalarType>& a_tangent_0,
    const NormalBase<ScalarType>& a_tangent_1, bool* a_requires_nudge,
    UnsignedIndex_t* a_split_counter, SurfaceOutputType* a_surface) {
  /* Defining constants and types */
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  using Plane = PlaneBase<ScalarType>;
  using AlignedParaboloid = AlignedParaboloidBase<ScalarType>;
  const ScalarType MACHINE_EPSILON = machine_epsilon<ScalarType>();
  const ScalarType DISTANCE_EPSILON = distance_epsilon<ScalarType>();
  const ScalarType ANGLE_EPSILON = angle_epsilon<ScalarType>();
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType SPLIT_TOL1 = static_cast<ScalarType>(0.9);
  const ScalarType SPLIT_TOL2 = static_cast<ScalarType>(0.999);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType FOUR = static_cast<ScalarType>(4);
  const ScalarType HALF = ONE / TWO;
  const ScalarType ONEQUARTER = ONE / FOUR;
  const ScalarType THREEQUARTERS = HALF + ONEQUARTER;
  const ScalarType ONE_HUNDRED = static_cast<ScalarType>(100);

  /* Function */
  const Pt& pt_ref = a_pt_ref.getPt();
  const Pt& pt_0 = a_pt_0.getPt();
  const Pt& pt_1 = a_pt_1.getPt();

  const Normal edge_vector = pt_1 - pt_0;
  Normal edge_vector_normalized = edge_vector;
  edge_vector_normalized.normalize();
  const ScalarType tgt0_dot_edge = a_tangent_0 * edge_vector_normalized;
  const ScalarType tgt1_dot_edge = a_tangent_1 * edge_vector_normalized;

  if (squaredMagnitude(edge_vector) < DISTANCE_EPSILON * DISTANCE_EPSILON &&
      fabs(ONE - tgt0_dot_edge) < ANGLE_EPSILON &&
      fabs(ONE + tgt1_dot_edge) < ANGLE_EPSILON) {
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      auto surface_arc = RationalBezierArc(
          pt_0.toDoublePt(), 0.5 * (pt_0.toDoublePt() + pt_1.toDoublePt()),
          pt_1.toDoublePt(), 0.0);
      surface_arc.reset_start_point_id(reinterpret_cast<std::uintptr_t>(&pt_0));
      surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_1));
      a_surface->addArc(surface_arc);
    }
    return ReturnType::fromScalarConstant(ZERO);
  }

  // bool split = ((a_tangent_0 * a_tangent_1) > SPLIT_TOL ||
  //               tgt0_dot_edge < ZERO || tgt1_dot_edge > ZERO);

  const ScalarType dot_tgt = a_tangent_0 * a_tangent_1;
  /* Split if tangents are both pointing towards a point far away */
  bool split = dot_tgt > SPLIT_TOL1;
  /* Split if the control point is on the wrong side of the edge */
  if (fabs(tgt0_dot_edge) < SPLIT_TOL2 && fabs(tgt1_dot_edge) < SPLIT_TOL2) {
    const Normal n_cross_t0 = crossProduct(a_plane_normal, a_tangent_0);
    const ScalarType triple_prod = n_cross_t0 * a_tangent_1;
    if (triple_prod * (n_cross_t0 * edge_vector) > ZERO) {
      split = true;
    }
    /* Else split if the tangents point opposite to the edge */
  } else if (tgt0_dot_edge < ZERO || tgt1_dot_edge > ZERO) {
    split = true;
  }

  if (split) {
    (*a_split_counter)++;
    const Pt average_pt = HALF * (pt_0 + pt_1);
    auto average_tangent = Normal(HALF * (a_tangent_0 + a_tangent_1));
    if (squaredMagnitude(average_tangent) < DISTANCE_EPSILON) {
      // average_tangent =
      //     Normal(ONEQUARTER * a_tangent_0 + THREEQUARTERS *
      //     a_tangent_1);
      average_tangent = crossProduct(edge_vector_normalized, a_plane_normal);
    }
    const ScalarType normal_correction = average_tangent * a_plane_normal;
    average_tangent = average_tangent - normal_correction * a_plane_normal;
    average_tangent.normalize();
    Pt projected_pt = projectPtAlongHalfLineOntoParaboloid<ScalarType>(
        a_paraboloid, average_tangent, average_pt);
    if (projected_pt[0] == static_cast<ScalarType>(DBL_MAX)) {
      *a_requires_nudge = true;
      return ReturnType::fromScalarConstant(ZERO);
    }
    if (*a_split_counter > 10) {
      if constexpr (std::is_same_v<ScalarType, double>) {
        *a_requires_nudge = true;
      }
      return ReturnType::fromScalarConstant(ZERO);
    }
    Normal tangent_projected_pt =
        computeAndCorrectTangentVectorAtPt<ScalarType>(
            a_paraboloid, a_plane_normal, pt_0, pt_1, a_tangent_1,
            projected_pt);

    // We need to store this vertex so that its address remains
    // unique over time
    if ((tangent_projected_pt[0] == ZERO && tangent_projected_pt[1] == ZERO &&
         tangent_projected_pt[2] == ZERO) ||
        squaredMagnitude(projected_pt - average_pt) <
            DISTANCE_EPSILON * DISTANCE_EPSILON) {
      if constexpr (!std::is_same_v<SurfaceOutputType, NoSurfaceOutput>) {
        auto surface_arc = RationalBezierArc(
            pt_0.toDoublePt(), 0.5 * (pt_0.toDoublePt() + pt_1.toDoublePt()),
            pt_1.toDoublePt(), 0.0);
        surface_arc.reset_start_point_id(
            reinterpret_cast<std::uintptr_t>(&pt_0));
        surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_1));
        a_surface->addArc(surface_arc);
      }
      return ReturnType::fromScalarConstant(ZERO);
    }

    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      Pt* new_point = new Pt(projected_pt);
      PtBase<double>* new_point_double =
          new PtBase<double>(static_cast<double>(projected_pt[0]),
                             static_cast<double>(projected_pt[1]),
                             static_cast<double>(projected_pt[2]));
      a_surface->addPt(new_point_double);
      return computeType3ContributionWithSplit<ReturnType, ScalarType>(
                 a_paraboloid, a_plane_normal, a_pt_ref, a_pt_0, *new_point,
                 a_tangent_0, tangent_projected_pt, a_requires_nudge,
                 a_split_counter, a_surface) +
             computeType3ContributionWithSplit<ReturnType, ScalarType>(
                 a_paraboloid, a_plane_normal, a_pt_ref, *new_point, a_pt_1,
                 -tangent_projected_pt, a_tangent_1, a_requires_nudge,
                 a_split_counter, a_surface);
    } else {
      // std::cout << "Split into " << pt_0 << " to " <<
      // projected_pt
      //           << " with tangents " << a_tangent_0 << " and "
      //           << tangent_projected_pt << std::endl;
      // std::cout << "Split into " << projected_pt << " to " <<
      // pt_1
      //           << " with tangents " << -tangent_projected_pt
      //           << " and "
      //           << a_tangent_1 << std::endl;

      return computeType3ContributionWithSplit<ReturnType, ScalarType>(
                 a_paraboloid, a_plane_normal, a_pt_ref, a_pt_0,
                 PtType(projected_pt), a_tangent_0, tangent_projected_pt,
                 a_requires_nudge, a_split_counter, a_surface) +
             computeType3ContributionWithSplit<ReturnType, ScalarType>(
                 a_paraboloid, a_plane_normal, a_pt_ref, PtType(projected_pt),
                 a_pt_1, -tangent_projected_pt, a_tangent_1, a_requires_nudge,
                 a_split_counter, a_surface);
    }
  } else {
    const auto arc = RationalBezierArcBase<ScalarType>(
        pt_0, a_tangent_0, pt_1, a_tangent_1, a_plane_normal, a_paraboloid);
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      auto surface_arc = RationalBezierArc(
          pt_0.toDoublePt(), a_tangent_0.toDoubleNormal(), pt_1.toDoublePt(),
          a_tangent_1.toDoubleNormal(), a_plane_normal.toDoubleNormal(),
          AlignedParaboloidBase<double>(a_paraboloid));
      surface_arc.reset_start_point_id(reinterpret_cast<std::uintptr_t>(&pt_0));
      surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_1));
      a_surface->addArc(surface_arc);
    }
    if (arc.weight() < ZERO) {
      *a_requires_nudge = true;
      return ReturnType::fromScalarConstant(ZERO);
    }
    auto moments =
        computeType3Contribution<ReturnType, ScalarType>(a_paraboloid, arc);
    if (!(&a_pt_ref == &a_pt_0 || &a_pt_ref == &a_pt_1)) {
      moments += computeTriangleCorrection<ReturnType, ScalarType>(
          a_paraboloid, pt_0, pt_1, pt_ref);
    }
    return moments;
  }
  // }
}

/************* Calculate arc contribution + gradient (with split)
 * *************/
template <class ReturnType, class ScalarType, class SurfaceOutputType,
          class PtType>
ReturnType computeType3ContributionWithGradientWithSplit(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtType& a_plane_normal, const PtType& a_pt_ref, const PtType& a_pt_0,
    const PtType& a_pt_1, const PtType& a_tangent_0, const PtType& a_tangent_1,
    SurfaceOutputType* a_surface) {
  /* Defining constants and types */
  using gradient_type = typename PtType::gradient_type;
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  using Plane = PlaneBase<ScalarType>;
  using AlignedParaboloid = AlignedParaboloidBase<ScalarType>;
  const ScalarType EPSILON =
      maximum(machine_epsilon<ScalarType>(), static_cast<ScalarType>(1.0e-24q));
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType FOUR = static_cast<ScalarType>(4);
  const ScalarType HALF = ONE / TWO;
  const ScalarType ONEQUARTER = ONE / FOUR;
  const ScalarType THREEQUARTERS = HALF + ONEQUARTER;
  const ScalarType ONE_HUNDRED = static_cast<ScalarType>(100);

  /* Function */
  const Pt& pt_ref = a_pt_ref.getPt();
  const Pt& pt_0 = a_pt_0.getPt();
  const Pt& pt_1 = a_pt_1.getPt();
  const Normal tgt_0 = Normal::fromPt(a_tangent_0.getPt());
  const Normal tgt_1 = Normal::fromPt(a_tangent_1.getPt());
  const Normal arc_normal = crossProduct(a_plane_normal.getPt(), pt_1 - pt_0);
  if ((tgt_0 * arc_normal > ZERO) != (tgt_1 * arc_normal > ZERO)) {
    std::cout << "Error: impossible tangent signs!" << std::endl;
    return ReturnType::fromScalarConstant(ZERO);
    // exit(-1);
  }
  if (squaredMagnitude(pt_1 - pt_0) < ONE_HUNDRED * EPSILON * EPSILON) {
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      auto surface_arc = RationalBezierArc(
          pt_1.toDoublePt(), 0.5 * (pt_0.toDoublePt() + pt_1.toDoublePt()),
          pt_0.toDoublePt(), 0.0);
      surface_arc.reset_start_point_id(reinterpret_cast<std::uintptr_t>(&pt_1));
      surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_0));
      a_surface->addArc(surface_arc);
    }
    return ReturnType::fromScalarConstant(ZERO);
  } else {
    const Normal edge_vector = pt_1 - pt_0;
    bool split = ((tgt_0 * edge_vector) < ZERO || (tgt_1 * edge_vector) > ZERO);
    if (split) {
      const auto average_pt_withgrad = HALF * (a_pt_0 + a_pt_1);
      auto avg_tgt_withgrad = HALF * (a_tangent_0 + a_tangent_1);
      if (squaredMagnitude(avg_tgt_withgrad.getPt()) <
          ONE_HUNDRED * EPSILON * EPSILON) {
        avg_tgt_withgrad =
            ONEQUARTER * a_tangent_0 + THREEQUARTERS * a_tangent_1;
      }
      Pt avg_tgt = avg_tgt_withgrad.getPt();
      const auto avg_tgt_grad = avg_tgt_withgrad.getData();
      const ScalarType norm_avg_tgt = magnitude(avg_tgt);
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
          computeAndCorrectTangentVectorAndGradientAtPt<ScalarType, PtType>(
              a_paraboloid, a_plane_normal, a_pt_0, a_pt_1, a_tangent_1,
              projected_pt_withgrad);
      // We need to store this vertex so that its address remains
      // unique over time
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        Pt* new_point = new Pt(projected_pt_withgrad.getPt());
        a_surface->addPt(new_point);
        auto new_proj_pt_withgrad = PtType(*new_point);
        new_proj_pt_withgrad.getData() = projected_pt_withgrad.getData();
        return computeType3ContributionWithGradientWithSplit<ReturnType,
                                                             ScalarType>(
                   a_paraboloid, a_plane_normal, a_pt_ref, a_pt_0,
                   new_proj_pt_withgrad, a_tangent_0,
                   tangent_projected_pt_withgrad, a_surface) +
               computeType3ContributionWithGradientWithSplit<ReturnType,
                                                             ScalarType>(
                   a_paraboloid, a_plane_normal, a_pt_ref, new_proj_pt_withgrad,
                   a_pt_1, -tangent_projected_pt_withgrad, a_tangent_1,
                   a_surface);
      } else {
        return computeType3ContributionWithGradientWithSplit<ReturnType,
                                                             ScalarType>(
                   a_paraboloid, a_plane_normal, a_pt_ref, a_pt_0,
                   projected_pt_withgrad, a_tangent_0,
                   tangent_projected_pt_withgrad, a_surface) +
               computeType3ContributionWithGradientWithSplit<ReturnType,
                                                             ScalarType>(
                   a_paraboloid, a_plane_normal, a_pt_ref,
                   projected_pt_withgrad, a_pt_1,
                   -tangent_projected_pt_withgrad, a_tangent_1, a_surface);
      }
    } else {
      const auto arc_with_gradient =
          RationalBezierArcWithGradientBase<PtType, ScalarType>(
              a_pt_0, a_tangent_0, a_pt_1, a_tangent_1, a_plane_normal,
              a_paraboloid);
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        auto surface_arc = RationalBezierArc(
            a_pt_0.getPt().toDoublePt(),
            arc_with_gradient.control_point().getPt().toDoublePt(),
            a_pt_1.getPt().toDoublePt(),
            static_cast<double>(arc_with_gradient.weight()));
        surface_arc.reset_start_point_id(
            reinterpret_cast<std::uintptr_t>(&(a_pt_0.getPt())));
        surface_arc.reset_end_point_id(
            reinterpret_cast<std::uintptr_t>(&(a_pt_1.getPt())));
        a_surface->addArc(surface_arc);
      }
      auto volume_with_gradient =
          computeType3ContributionWithGradient<ReturnType, ScalarType>(
              a_paraboloid, arc_with_gradient);
      if (!(&a_pt_ref == &a_pt_0 || &a_pt_ref == &a_pt_1)) {
        volume_with_gradient +=
            computeTriangleCorrectionWithGradient<ReturnType, ScalarType>(
                a_paraboloid, a_pt_0, a_pt_1, a_pt_ref);
      }
      return volume_with_gradient;
    }
  }
}

/**************** Calculate line contribution ******************/
// Starts from an entry, returns the exit that is reached.
template <class ReturnType, class ScalarType, class HalfEdgeType, class PtType>
ReturnType computeUnclippedSegmentType1Contribution(
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid,
    const PtType& a_ref_pt, const HalfEdgeType a_entry_half_edge,
    HalfEdgeType& a_exit_half_edge, const bool skip_first) {
  /* Defining constants and types */
  const ScalarType ZERO = static_cast<ScalarType>(0);

  /* Function */
  ReturnType full_moments = ReturnType::fromScalarConstant(ZERO);

  assert(a_entry_half_edge->getPreviousVertex()->isClipped() ||
         a_entry_half_edge->getPreviousVertex()->doesNotNeedToSeek());
  assert(a_entry_half_edge->getNextHalfEdge()->getVertex()->isNotClipped());

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
      if constexpr (!has_embedded_gradient<ReturnType>::value) {
        full_moments += computeType1Contribution<ReturnType, ScalarType>(
            a_ref_pt, prev_pt, curr_pt);
      } else {
        full_moments +=
            computeType1ContributionWithGradient<ReturnType, ScalarType>(
                a_ref_pt, prev_pt, curr_pt);
      }
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

/**************** Calculate all arc contributions ******************/
template <class ReturnType, class ScalarType, class HalfEdgeType,
          class SurfaceOutputType, class PtType>
ReturnType computeNewEdgeSegmentContribution(
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid,
    const PtType& a_ref_pt, const HalfEdgeType a_entry_half_edge,
    const HalfEdgeType a_exit_half_edge, const bool skip_first,
    const bool a_ignore_type3, bool* a_requires_nudge,
    SurfaceOutputType* a_surface) {
  ReturnType full_moments =
      ReturnType::fromScalarConstant(static_cast<ScalarType>(0.0));
  // Handle new edge on exit->entry
  if constexpr (!has_embedded_gradient<ReturnType>::value) {
    if (!skip_first) {
      full_moments += computeType1Contribution<ReturnType, ScalarType>(
          a_ref_pt, a_exit_half_edge->getVertex()->getLocation(),
          a_entry_half_edge->getVertex()->getLocation());
    }
    full_moments += computeType2Contribution<ReturnType, ScalarType>(
        a_aligned_paraboloid, a_exit_half_edge->getVertex()->getLocation(),
        a_entry_half_edge->getVertex()->getLocation());
  } else {
    if (!skip_first) {
      full_moments +=
          computeType1ContributionWithGradient<ReturnType, ScalarType>(
              a_ref_pt, a_exit_half_edge->getVertex()->getLocation(),
              a_entry_half_edge->getVertex()->getLocation());
    }
    full_moments +=
        computeType2ContributionWithGradient<ReturnType, ScalarType>(
            a_aligned_paraboloid, a_exit_half_edge->getVertex()->getLocation(),
            a_entry_half_edge->getVertex()->getLocation());
  }

  if (!a_ignore_type3) {
    if constexpr (has_embedded_gradient<ReturnType>::value) {
      full_moments +=
          orientAndApplyType3CorrectionWithGradients<ReturnType, ScalarType>(
              a_aligned_paraboloid, a_exit_half_edge, a_entry_half_edge,
              a_surface, a_requires_nudge);
    } else {
      full_moments += orientAndApplyType3Correction<ReturnType, ScalarType>(
          a_aligned_paraboloid, a_exit_half_edge, a_entry_half_edge,
          a_requires_nudge, a_surface);
    }
  }
  return full_moments;
}

/************* Calculate moments from non-aligned paraboloid **********/
template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class ParaboloidType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithParaboloid(SegmentedHalfEdgePolyhedronType* a_polytope,
                                  HalfEdgePolytopeType* a_complete_polytope,
                                  const ParaboloidType& a_paraboloid) {
  /* Defining type aliases (needed to ensure precision is
   * consistent)*/
  using ScalarType = typename ParaboloidType::value_type;
  using PtType = typename SegmentedHalfEdgePolyhedronType::pt_type;
  static_assert(std::is_same_v<typename PtType::value_type, ScalarType>);
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  using Plane = PlaneBase<ScalarType>;
  using AlignedParaboloid = AlignedParaboloidBase<ScalarType>;

  /* Defining constants */
  const ScalarType DISTANCE_EPSILON = distance_epsilon<ScalarType>();
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType THREE = static_cast<ScalarType>(3);

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
      moments.getMoments() = ReturnType::moment_type::fromScalarConstant(ZERO);
      moments.getSurface().setParaboloid(a_paraboloid);
    } else {
      moments = ReturnType::fromScalarConstant(ZERO);
    }
    return moments;
  }

  // Move into reconstruction reference frame
  const UnsignedIndex_t original_number_of_vertices =
      a_polytope->getNumberOfVertices();
  const auto& datum = a_paraboloid.getDatum();
  const auto& ref_frame = a_paraboloid.getReferenceFrame();
  // assert(ref_frame.isOrthonormalBasis());

  const Pt start_pt = a_polytope->getVertex(0)->getLocation().getPt() - datum;
  ScalarType max_dist_sq = ZERO;
  for (UnsignedIndex_t v = 0; v < original_number_of_vertices; ++v) {
    const Pt original_pt =
        a_polytope->getVertex(v)->getLocation().getPt() - datum;
    if (v > 0) {
      max_dist_sq =
          maximum(max_dist_sq, squaredMagnitude(original_pt - start_pt));
    }
    PtType projected_location;
    auto& pt = projected_location.getPt();
    for (UnsignedIndex_t n = 0; n < 3; ++n) {
      pt[n] = ref_frame[n] * original_pt;
    }
    a_polytope->getVertex(v)->setLocation(projected_location);
  }

  const ScalarType inv_scale = maximum(
      static_cast<ScalarType>(1.0e6) * DISTANCE_EPSILON, sqrt(max_dist_sq));
  const ScalarType inv_volume_scale = inv_scale * max_dist_sq;
  const ScalarType scale = ONE / inv_scale;
  const ScalarType volume_scale = scale * scale * scale;

  // Normalized polyhedron
  for (UnsignedIndex_t v = 0; v < original_number_of_vertices; ++v) {
    auto& pt = a_polytope->getVertex(v)->getLocation().getPt();
    pt *= scale;
  }
  // Normalized paraboloid
  auto const& scaled_aligned_paraboloid =
      AlignedParaboloid(std::array<ScalarType, 2>{
          a_paraboloid.getAlignedParaboloid().a() * inv_scale,
          a_paraboloid.getAlignedParaboloid().b() * inv_scale});

  // // Recalculate in face plane information
  // for (auto& face : (*a_polytope)) {
  //   auto normal = Normal(ZERO, ZERO, ZERO);
  //   const auto starting_half_edge =
  //   face->getStartingHalfEdge(); auto current_half_edge =
  //   starting_half_edge; auto next_half_edge =
  //   starting_half_edge->getNextHalfEdge(); const auto&
  //   start_location =
  //       starting_half_edge->getPreviousVertex()->getLocation().getPt();
  //   do {
  //     normal += crossProduct(
  //         current_half_edge->getVertex()->getLocation().getPt()
  //         -
  //             start_location,
  //         next_half_edge->getVertex()->getLocation().getPt() -
  //         start_location);
  //     current_half_edge = next_half_edge;
  //     next_half_edge = next_half_edge->getNextHalfEdge();
  //   } while (next_half_edge != starting_half_edge);
  //   normal.normalize();
  //   if (squaredMagnitude(normal) <
  //   static_cast<ScalarType>(0.99)) {
  //     normal = Normal(ONE / THREE, ONE / THREE, ONE / THREE);
  //   }
  //   face->setPlane(Plane(normal, normal * start_location));
  // }

  // Compute intersection
  if constexpr (has_paraboloid_surface<ReturnType>::value) {
    moments.getSurface().setParaboloid(a_paraboloid);
    moments.getMoments() = intersectPolyhedronWithAlignedParaboloid<
        typename ReturnType::moment_type>(a_polytope, a_complete_polytope,
                                          scaled_aligned_paraboloid,
                                          &moments.getSurface());
  } else {
    NoSurfaceOutput* surf = nullptr;
    moments = intersectPolyhedronWithAlignedParaboloid<ReturnType>(
        a_polytope, a_complete_polytope, scaled_aligned_paraboloid, surf);
  }

  // bool error = false;
  // if constexpr (!has_paraboloid_surface<ReturnType>::value &&
  //               !has_embedded_gradient<ReturnType>::value) {
  //   if constexpr ((std::is_same_v<ReturnType, VolumeMomentsBase<double>>
  //   ||
  //                  std::is_same_v<ReturnType,
  //                  VolumeMomentsBase<Quad_t>>)) {
  //     if (moments.volume() < -999.9) {
  //       error = true;
  //     }
  //   } else {
  //     if (moments < -999.9) {  // || fabs(moments.volume() *
  //                              // inv_volume_scale -
  //                              //     0.0036361026083212) < 1.0e-9) {
  //       error = true;
  //     }
  //   }
  // }

  // if (error) {
  //   sleep(1);
  //   std::cout << "SCALE = " << inv_scale << std::endl;
  //   std::cout << "PARABOLOID = " << scaled_aligned_paraboloid <<
  //   std::endl; std::cout << (*a_polytope) << std::endl; const
  //   UnsignedIndex_t number_of_vertices =
  //       a_polytope->getNumberOfVertices();
  //   for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
  //     auto& pt = a_polytope->getVertex(v)->getLocation().getPt();
  //     pt *= inv_scale;
  //   }
  //   std::cout << (*a_polytope) << std::endl;
  //   // exit(-1);
  // }

  // Un-normalized polyhedron
  // const UnsignedIndex_t number_of_vertices =
  // a_polytope->getNumberOfVertices(); for (UnsignedIndex_t v = 0; v <
  // number_of_vertices; ++v) {
  //   auto& pt = a_polytope->getVertex(v)->getLocation().getPt();
  //   pt *= inv_scale;
  // }

  // Un-normalized moments
  if constexpr (has_paraboloid_surface<ReturnType>::value) {
    if constexpr ((std::is_same_v<typename ReturnType::moment_type,
                                  VolumeMomentsBase<double>> ||
                   std::is_same_v<typename ReturnType::moment_type,
                                  VolumeMomentsBase<Quad_t>>)) {
      moments.getMoments().centroid().getPt() *= inv_volume_scale * inv_scale;
      moments.getMoments().volume() *= inv_volume_scale;
    } else {
      moments.getMoments().volume() *= inv_volume_scale;
    }
    auto& arc_list = moments.getSurface().getArcs();
    for (std::size_t i = 0; i < arc_list.size(); ++i) {
      arc_list[i].start_point() *= static_cast<double>(inv_scale);
      arc_list[i].control_point() *= static_cast<double>(inv_scale);
      arc_list[i].end_point() *= static_cast<double>(inv_scale);
    }
  } else if constexpr ((std::is_same_v<ReturnType, VolumeMomentsBase<double>> ||
                        std::is_same_v<ReturnType,
                                       VolumeMomentsBase<Quad_t>>)) {
    moments.centroid().getPt() *= inv_volume_scale * inv_scale;
    moments.volume() *= inv_volume_scale;
  } else {
    moments.volume() *= inv_volume_scale;
  }

  // Un-normalized gradient
  if constexpr (has_embedded_gradient<ReturnType>::value) {
    if constexpr ((std::is_same_v<ReturnType, VolumeMomentsBase<double>> ||
                   std::is_same_v<ReturnType, VolumeMomentsBase<Quad_t>>)) {
      const auto normed_grad_m1 = moments.centroid().getData();
      const auto normed_grad_m0 = moments.volume_gradient();
      auto& grad_m1 = moments.centroid().getData();
      auto& grad_m0 = moments.volume_gradient();
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        grad_m1[d].setGradRx(normed_grad_m1[d].getGradRx() * inv_volume_scale *
                             inv_scale);
        grad_m1[d].setGradRy(normed_grad_m1[d].getGradRy() * inv_volume_scale *
                             inv_scale);
        grad_m1[d].setGradRz(normed_grad_m1[d].getGradRz() * inv_volume_scale *
                             inv_scale);
        grad_m1[d].setGradTx(normed_grad_m1[d].getGradTx() * inv_volume_scale);
        grad_m1[d].setGradTy(normed_grad_m1[d].getGradTy() * inv_volume_scale);
        grad_m1[d].setGradTz(normed_grad_m1[d].getGradTz() * inv_volume_scale);
        grad_m1[d].setGradA(normed_grad_m1[d].getGradA() * inv_volume_scale *
                            inv_scale * inv_scale);
        grad_m1[d].setGradB(normed_grad_m1[d].getGradB() * inv_volume_scale *
                            inv_scale * inv_scale);
      }
      grad_m0.setGradRx(normed_grad_m0.getGradRx() * inv_volume_scale);
      grad_m0.setGradRy(normed_grad_m0.getGradRy() * inv_volume_scale);
      grad_m0.setGradRz(normed_grad_m0.getGradRz() * inv_volume_scale);
      grad_m0.setGradTx(normed_grad_m0.getGradTx() * inv_scale * inv_scale);
      grad_m0.setGradTy(normed_grad_m0.getGradTy() * inv_scale * inv_scale);
      grad_m0.setGradTz(normed_grad_m0.getGradTz() * inv_scale * inv_scale);
      grad_m0.setGradA(normed_grad_m0.getGradA() * inv_volume_scale *
                       inv_scale);
      grad_m0.setGradB(normed_grad_m0.getGradB() * inv_volume_scale *
                       inv_scale);
    } else {
      const auto normed_grad_m0 = moments.volume_gradient();
      auto& grad_m0 = moments.volume_gradient();
      grad_m0.setGradRx(normed_grad_m0.getGradRx() * inv_volume_scale);
      grad_m0.setGradRy(normed_grad_m0.getGradRy() * inv_volume_scale);
      grad_m0.setGradRz(normed_grad_m0.getGradRz() * inv_volume_scale);
      grad_m0.setGradTx(normed_grad_m0.getGradTx() * inv_scale * inv_scale);
      grad_m0.setGradTy(normed_grad_m0.getGradTy() * inv_scale * inv_scale);
      grad_m0.setGradTz(normed_grad_m0.getGradTz() * inv_scale * inv_scale);
      grad_m0.setGradA(normed_grad_m0.getGradA() * inv_volume_scale *
                       inv_scale);
      grad_m0.setGradB(normed_grad_m0.getGradB() * inv_volume_scale *
                       inv_scale);
    }
  }

  // Rotate base polyhedron back
  // for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
  //   const Pt& original_pt =
  //   a_polytope->getVertex(v)->getLocation().getPt(); typename
  //   SegmentedHalfEdgePolyhedronType::pt_type projected_location; auto& pt
  //   = projected_location.getPt(); pt = Pt(ZERO, ZERO, ZERO); for
  //   (UnsignedIndex_t d = 0; d < 3; ++d) {
  //     for (UnsignedIndex_t n = 0; n < 3; ++n) {
  //       pt[n] += ref_frame[d][n] * original_pt[d];
  //     }
  //   }
  //   pt += datum;
  //   a_polytope->getVertex(v)->setLocation(projected_location);
  // }

  // Rotate gradient of first moment
  if constexpr (has_embedded_gradient<ReturnType>::value) {
    if constexpr ((std::is_same_v<typename ReturnType::moment_type,
                                  VolumeMomentsBase<double>> ||
                   std::is_same_v<typename ReturnType::moment_type,
                                  VolumeMomentsBase<Quad_t>>)) {
      using gradient_type = typename ReturnType::gradient_type;
      // Compute gradient of ref frame
      std::array<std::array<gradient_type, 3>, 3> frame_grad;
      std::array<gradient_type, 3> datum_grad;
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        datum_grad[d] = gradient_type(ZERO);
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          frame_grad[d][n] = gradient_type(ZERO);
        }
      }
      frame_grad[0][0].setGradRy(ref_frame[0][2] * ref_frame[1][1] -
                                 ref_frame[0][1] * ref_frame[1][2]);
      frame_grad[0][1].setGradRy(-ref_frame[0][2] * ref_frame[1][0] +
                                 ref_frame[0][0] * ref_frame[1][2]);
      frame_grad[0][2].setGradRy(ref_frame[0][1] * ref_frame[1][0] -
                                 ref_frame[0][0] * ref_frame[1][1]);
      frame_grad[0][0].setGradRz(ref_frame[0][2] * ref_frame[2][1] -
                                 ref_frame[0][1] * ref_frame[2][2]);
      frame_grad[0][1].setGradRz(-ref_frame[0][2] * ref_frame[2][0] +
                                 ref_frame[0][0] * ref_frame[2][2]);
      frame_grad[0][2].setGradRz(ref_frame[0][1] * ref_frame[2][0] -
                                 ref_frame[0][0] * ref_frame[2][1]);
      frame_grad[1][0].setGradRx(-ref_frame[0][2] * ref_frame[1][1] +
                                 ref_frame[0][1] * ref_frame[1][2]);
      frame_grad[1][1].setGradRx(ref_frame[0][2] * ref_frame[1][0] -
                                 ref_frame[0][0] * ref_frame[1][2]);
      frame_grad[1][2].setGradRx(-ref_frame[0][1] * ref_frame[1][0] +
                                 ref_frame[0][0] * ref_frame[1][1]);
      frame_grad[1][0].setGradRz(ref_frame[1][2] * ref_frame[2][1] -
                                 ref_frame[1][1] * ref_frame[2][2]);
      frame_grad[1][1].setGradRz(-ref_frame[1][2] * ref_frame[2][0] +
                                 ref_frame[1][0] * ref_frame[2][2]);
      frame_grad[1][2].setGradRz(ref_frame[1][1] * ref_frame[2][0] -
                                 ref_frame[1][0] * ref_frame[2][1]);
      frame_grad[2][0].setGradRx(-ref_frame[0][2] * ref_frame[2][1] +
                                 ref_frame[0][1] * ref_frame[2][2]);
      frame_grad[2][1].setGradRx(ref_frame[0][2] * ref_frame[2][0] -
                                 ref_frame[0][0] * ref_frame[2][2]);
      frame_grad[2][2].setGradRx(-ref_frame[0][1] * ref_frame[2][0] +
                                 ref_frame[0][0] * ref_frame[2][1]);
      frame_grad[2][0].setGradRy(-ref_frame[1][2] * ref_frame[2][1] +
                                 ref_frame[1][1] * ref_frame[2][2]);
      frame_grad[2][1].setGradRy(ref_frame[1][2] * ref_frame[2][0] -
                                 ref_frame[1][0] * ref_frame[2][2]);
      frame_grad[2][2].setGradRy(-ref_frame[1][1] * ref_frame[2][0] +
                                 ref_frame[1][0] * ref_frame[2][1]);
      datum_grad[0].setGradTx(ref_frame[0][0]);
      datum_grad[0].setGradTy(ref_frame[1][0]);
      datum_grad[0].setGradTz(ref_frame[2][0]);
      datum_grad[1].setGradTx(ref_frame[0][1]);
      datum_grad[1].setGradTy(ref_frame[1][1]);
      datum_grad[1].setGradTz(ref_frame[2][1]);
      datum_grad[2].setGradTx(ref_frame[0][2]);
      datum_grad[2].setGradTy(ref_frame[1][2]);
      datum_grad[2].setGradTz(ref_frame[2][2]);
      auto pt_grad = std::array<gradient_type, 3>(
          {gradient_type(ZERO), gradient_type(ZERO), gradient_type(ZERO)});
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          pt_grad[n] += ref_frame[d][n] * moments.centroid().getData()[d] +
                        frame_grad[d][n] * moments.centroid().getPt()[d];
        }
      }
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        pt_grad[d] += moments.volume_gradient() * datum[d] +
                      moments.volume() * datum_grad[d];
      }
      moments.centroid().getData() = pt_grad;

      // Move first moment back to original frame of reference
      auto pt = Pt(ZERO, ZERO, ZERO);
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          pt[n] += ref_frame[d][n] * moments.centroid().getPt()[d];
        }
      }
      pt += moments.volume() * datum;
      moments.centroid().getPt() = pt;
    }
  }
  // Move first moment back to original frame of reference
  if constexpr (has_paraboloid_surface<ReturnType>::value) {
    if constexpr (std::is_same_v<typename ReturnType::moment_type,
                                 VolumeMomentsBase<double>> ||
                  std::is_same_v<typename ReturnType::moment_type,
                                 VolumeMomentsBase<Quad_t>>) {
      auto pt = Pt(ZERO, ZERO, ZERO);
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          pt[n] += ref_frame[d][n] * moments.getMoments().centroid().getPt()[d];
        }
      }
      pt += moments.getMoments().volume() * datum;
      moments.getMoments().centroid().getPt() = pt;
    }
  }
  if constexpr (!has_paraboloid_surface<ReturnType>::value &&
                (std::is_same_v<ReturnType, VolumeMomentsBase<double>> ||
                 std::is_same_v<ReturnType, VolumeMomentsBase<Quad_t>>)) {
    auto pt = Pt(ZERO, ZERO, ZERO);
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      for (UnsignedIndex_t n = 0; n < 3; ++n) {
        pt[n] += ref_frame[d][n] * moments.centroid().getPt()[d];
      }
    }
    pt += moments.volume() * datum;
    moments.centroid().getPt() = pt;
  }

  return moments;
}

/*********** Calculate moments from aligned paraboloid *************/
template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class AlignedParaboloidType,
          class SurfaceOutputType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithAlignedParaboloid(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AlignedParaboloidType& a_paraboloid, SurfaceOutputType* a_surface) {
  // Assumed a_polytope is already rotated to be in same
  // coordinate system as a_paraboloid.

  // Below function computes the entire integration
  const auto moments = formParaboloidIntersectionBases<ReturnType>(
      a_polytope, a_complete_polytope, a_paraboloid, 0, a_surface);

  return moments;
}

/******************* Place one intersection on segment **********************/
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

/******************* Place two intersection on segment **********************/
template <class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType, class VertexType>
void placeDoubleIntercept(
    const StackVector<VertexType, 2>& a_intersection_location,
    HalfEdgeType* a_half_edge_with_intersection,
    SegmentedHalfEdgePolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope) {
  assert(a_intersection_location.size() == 2);

  // Need to place furthest vertex first so that
  // a_half_edge_with_intersection remains attached to same vertex
  // it started on. Also want to keep property that the half edge
  // the vertex stores has the previousVertex also unclipped. This
  // means for the first (furthest) vertex, need to reference the
  // opposite half edge of the current one. This enables more
  // efficient face truncation in the next phase of the algorithm,
  // since new vertices will always have a half edge going from
  // unclipped->new->clipped.

  placeSingleIntercept(a_intersection_location[0],
                       a_half_edge_with_intersection, a_polytope,
                       a_complete_polytope);
  a_half_edge_with_intersection->getPreviousVertex()->setHalfEdge(
      a_half_edge_with_intersection->getOppositeHalfEdge());
  placeSingleIntercept(a_intersection_location[1],
                       a_half_edge_with_intersection, a_polytope,
                       a_complete_polytope);
}

/******************* Find intersections on segment **********************/
template <class PtType, class ScalarType>
void checkAndFindIntercepts(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid, const PtType& a_pt_0,
    const PtType& a_pt_1, StackVector<PtType, 2>* a_intercepts,
    const ScalarType a_nudge_epsilon, const bool a_elliptic);

template <class PtType, class ScalarType>
inline enable_if_t<!has_embedded_gradient<PtType>::value, void>
checkAndFindIntercepts(const AlignedParaboloidBase<ScalarType>& a_paraboloid,
                       const PtBase<ScalarType>& a_pt_0,
                       const PtBase<ScalarType>& a_pt_1,
                       StackVector<PtBase<ScalarType>, 2>* a_intercepts,
                       const ScalarType a_nudge_epsilon,
                       const bool a_elliptic) {
  static_assert(std::is_same_v<PtType, PtBase<ScalarType>>);
  const ScalarType EPSILON_LO = -static_cast<ScalarType>(0.5) * a_nudge_epsilon;
  const ScalarType EPSILON_HI = static_cast<ScalarType>(1) - EPSILON_LO;
  a_intercepts->resize(0);
  // bool solve = true;

  // /* Shortcut for elliptic case */
  // if (a_elliptic) {
  //   if (a_paraboloid.a() > ZERO) {
  //     if (vertexBelow(a_pt_0, a_paraboloid) &&
  //         vertexBelow(a_pt_1, a_paraboloid)) {
  //       solve = false;
  //     }
  //   } else if (!vertexBelow(a_pt_0, a_paraboloid) &&
  //              !vertexBelow(a_pt_1, a_paraboloid)) {
  //     solve = false;
  //   }
  // }

  // /* Compute potential intersection */
  // if (solve) {
  const auto& pt_0 = a_pt_0.getPt();
  const auto& pt_1 = a_pt_1.getPt();
  const auto pt_diff = pt_1 - pt_0;
  const ScalarType a = a_paraboloid.a() * pt_diff[0] * pt_diff[0] +
                       a_paraboloid.b() * pt_diff[1] * pt_diff[1];
  const ScalarType b =
      static_cast<ScalarType>(2) * (a_paraboloid.a() * pt_diff[0] * pt_0[0] +
                                    a_paraboloid.b() * pt_diff[1] * pt_0[1]) +
      pt_diff[2];
  const ScalarType c = a_paraboloid.a() * pt_0[0] * pt_0[0] +
                       a_paraboloid.b() * pt_0[1] * pt_0[1] + pt_0[2];

  const StackVector<ScalarType, 2> solutions =
      solveQuadratic<ScalarType>(a, b, c);

  for (auto& solution : solutions) {
    if (solution > EPSILON_LO && solution < EPSILON_HI) {
      // if (solution > ZERO && solution < ONE) {
      a_intercepts->push_back(PtBase<ScalarType>(pt_0 + solution * pt_diff));
    }
  }

  // if (solutions.size() > 0) {
  //   if (solutions[0] >= ZERO) {
  //     if (solutions[1] <= ONE) {
  //       a_intercepts->push_back(
  //           PtBase<ScalarType>(pt_0 + solutions[0] *
  //           pt_diff));
  //       a_intercepts->push_back(
  //           PtBase<ScalarType>(pt_0 + solutions[1] *
  //           pt_diff));
  //     } else if (solutions[0] <= ONE) {
  //       a_intercepts->push_back(
  //           PtBase<ScalarType>(pt_0 + solutions[0] *
  //           pt_diff));
  //     }
  //   } else if (solutions[1] >= ZERO && solutions[1] <= ONE) {
  //     a_intercepts->push_back(
  //         PtBase<ScalarType>(pt_0 + solutions[1] * pt_diff));
  //   }
  // }
  // }
}

template <class PtType, class ScalarType>
inline enable_if_t<has_embedded_gradient<PtType>::value, void>
checkAndFindIntercepts(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtWithGradient<typename PtType::gradient_type>& a_pt_0,
    const PtWithGradient<typename PtType::gradient_type>& a_pt_1,
    StackVector<PtWithGradient<typename PtType::gradient_type>, 2>*
        a_intercepts,
    const ScalarType a_nudge_epsilon, const bool a_elliptic) {
  /* Defining constants and types */
  using gradient_type = typename PtType::gradient_type;
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);

  /* Function */
  const auto pt_diff_with_grad = a_pt_1 - a_pt_0;
  const auto& pt_0 = a_pt_0.getPt();
  const auto& pt_1 = a_pt_1.getPt();
  const auto& pt_diff = pt_diff_with_grad.getPt();
  const auto& pt_0_grad = a_pt_0.getData();
  const auto& pt_diff_grad = pt_diff_with_grad.getData();
  const ScalarType A = a_paraboloid.a(), B = a_paraboloid.b();
  auto A_grad = gradient_type(ZERO), B_grad = gradient_type(ZERO);
  A_grad.setGradA(ONE);
  B_grad.setGradB(ONE);
  const ScalarType a =
      A * pt_diff[0] * pt_diff[0] + B * pt_diff[1] * pt_diff[1];
  const ScalarType b =
      TWO * (A * pt_diff[0] * pt_0[0] + B * pt_diff[1] * pt_0[1]) + pt_diff[2];
  const ScalarType c = A * pt_0[0] * pt_0[0] + B * pt_0[1] * pt_0[1] + pt_0[2];
  const auto a_grad = A_grad * pt_diff[0] * pt_diff[0] +
                      B_grad * pt_diff[1] * pt_diff[1] +
                      A * (TWO * pt_diff_grad[0] * pt_diff[0]) +
                      B * (TWO * pt_diff_grad[1] * pt_diff[1]);
  const auto b_grad =
      TWO * (A_grad * pt_diff[0] * pt_0[0] + B_grad * pt_diff[1] * pt_0[1] +
             A * (pt_diff_grad[0] * pt_0[0] + pt_diff[0] * pt_0_grad[0]) +
             B * (pt_diff_grad[1] * pt_0[1] + pt_diff[1] * pt_0_grad[1])) +
      pt_diff_grad[2];
  const auto c_grad = A_grad * pt_0[0] * pt_0[0] + B_grad * pt_0[1] * pt_0[1] +
                      A * TWO * pt_0_grad[0] * pt_0[0] +
                      B * TWO * pt_0_grad[1] * pt_0[1] + pt_0_grad[2];
  a_intercepts->resize(0);
  const auto solutions = solveQuadraticWithGradient<gradient_type>(
      a, b, c, a_grad, b_grad, c_grad);
  for (const auto& solution : solutions) {
    const auto sol = solution.first;
    const auto sol_grad = solution.second;
    if (sol >= -a_nudge_epsilon && sol <= ONE + a_nudge_epsilon) {
      auto intersection = PtType(Pt(pt_0 + sol * pt_diff));
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        intersection.getData()[d] =
            pt_0_grad[d] + sol_grad * pt_diff[d] + sol * pt_diff_grad[d];
      }
      a_intercepts->push_back(PtType(intersection));
    }
  }
}

/**************** Flag: is vertex below aligned paraboloid?
 * *****************/
template <class VertexType>
bool vertexBelow(const VertexType& a_pt,
                 const AlignedParaboloidBase<typename VertexType::value_type>&
                     a_paraboloid) {
  const auto& pt = a_pt.getPt();
  return pt[2] <
         -(a_paraboloid.a() * pt[0] * pt[0] + a_paraboloid.b() * pt[1] * pt[1]);
}

// This algorithm is based on the method shared from the website below, and
// from several stackoverflow posts citing that website.
// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
template <class ScalarType>
inline bool isPtBeforeIntersectionWithEdge(
    const std::array<ScalarType, 2>& a_test_pt,
    const PtBase<ScalarType>& a_vertex_0,
    const PtBase<ScalarType>& a_vertex_1) {
  if ((a_test_pt[1] > a_vertex_0[1]) == (a_test_pt[1] > a_vertex_1[1])) {
    return false;  // Projected ray never intersects edge.
  }
  const ScalarType location_of_intersection_along_ray =
      (a_vertex_0[0] - a_vertex_1[0]) * (a_test_pt[1] - a_vertex_1[1]) /
          (a_vertex_0[1] - a_vertex_1[1]) +
      a_vertex_1[0];
  // Intersection was to the right if
  // location_of_intersection_along_ray is greater.
  return a_test_pt[0] < location_of_intersection_along_ray;
}

// This algorithm is based on the method shared from the website below, and
// from several stackoverflow posts citing that website.
// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
template <class ScalarType>
inline bool isPtBeforeIntersectionWithEdgeWithComponent(
    const PtBase<ScalarType>& a_test_pt, const PtBase<ScalarType>& a_vertex_0,
    const PtBase<ScalarType>& a_vertex_1, const UnsignedIndex_t a_index) {
  const UnsignedIndex_t id0 = a_index;
  const UnsignedIndex_t id1 = (a_index + 1) % 3;
  if ((a_test_pt[id1] > a_vertex_0[id1]) ==
      (a_test_pt[id1] > a_vertex_1[id1])) {
    return false;  // Projected ray never intersects edge.
  }
  const ScalarType location_of_intersection_along_ray =
      (a_vertex_0[id0] - a_vertex_1[id0]) * (a_test_pt[id1] - a_vertex_1[id1]) /
          (a_vertex_0[id1] - a_vertex_1[id1]) +
      a_vertex_1[id0];
  // Intersection was to the right if
  // location_of_intersection_along_ray is greater.
  return a_test_pt[id0] < location_of_intersection_along_ray;
}

// If centroid is outside of polygon, or any vertex on face
// is inside the ellipse, then the ellipse is not contained by the face.
template <class ScalarType, class HalfEdgeType>
bool ellipseContainedInFace(
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid,
    const PlaneBase<ScalarType>& a_face_plane,
    HalfEdgeType* const a_half_edge) {
  /* Defining constants and types */
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType TWO = static_cast<ScalarType>(2);

  /* Function */
  const auto& face_normal = a_face_plane.normal();
  const std::array<ScalarType, 2> conic_center{
      {face_normal[0] / (TWO * a_aligned_paraboloid.a() * face_normal[2]),
       face_normal[1] / (TWO * a_aligned_paraboloid.b() * face_normal[2])}};
  const ScalarType delta_face = a_face_plane.distance() / face_normal[2];
  const ScalarType gamma_face =
      a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
      a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] - delta_face;
  if (a_aligned_paraboloid.a() * gamma_face < ZERO) {
    return false;
  }

  // First we will check if centroid is in the bounding box
  // of the face polygon. Due to the fact we are assuming
  // the paraboloid was in Z direction, we assume the ellipse
  // lives on the x/y plane and project the polygon
  // down to it as well (essentially neglecting the z component).
  auto current_half_edge = a_half_edge;
  std::array<ScalarType, 2> xy_min{
      {static_cast<ScalarType>(DBL_MAX), static_cast<ScalarType>(DBL_MAX)}};
  std::array<ScalarType, 2> xy_max{
      {-static_cast<ScalarType>(DBL_MAX), -static_cast<ScalarType>(DBL_MAX)}};
  do {
    const PtBase<ScalarType>& location =
        current_half_edge->getVertex()->getLocation().getPt();
    for (UnsignedIndex_t d = 0; d < 2; ++d) {
      xy_min[d] = minimum(xy_min[d], location[d]);
      xy_max[d] = maximum(xy_max[d], location[d]);
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
    const PtBase<ScalarType>& location_0 =
        current_half_edge->getPreviousVertex()->getLocation().getPt();
    const PtBase<ScalarType>& location_1 =
        current_half_edge->getVertex()->getLocation().getPt();
    if (isPtBeforeIntersectionWithEdge<ScalarType>(conic_center, location_0,
                                                   location_1)) {
      pt_internal_to_polygon = !pt_internal_to_polygon;
    }
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != a_half_edge);
  return pt_internal_to_polygon;
}

template <class ReturnType, class ScalarType, class HalfEdgeType,
          class SurfaceOutputType>
ReturnType orientAndApplyType3Correction(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    HalfEdgeType* a_start, HalfEdgeType* a_end, bool* a_requires_nudge,
    SurfaceOutputType* a_surface) {
  /* Defining constants and types */
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  using Plane = PlaneBase<ScalarType>;
  using AlignedParaboloid = AlignedParaboloidBase<ScalarType>;
  const ScalarType MACHINE_EPSILON = machine_epsilon<ScalarType>();
  const ScalarType DISTANCE_EPSILON = distance_epsilon<ScalarType>();
  const ScalarType ANGLE_EPSILON = angle_epsilon<ScalarType>();
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType HALF = ONE / TWO;
  const ScalarType TEN = static_cast<ScalarType>(10);
  const ScalarType ONE_HUNDRED = static_cast<ScalarType>(100);

  /* Function */
  const auto& pt_0 = a_start->getVertex()->getLocation().getPt();
  const auto& pt_1 = a_end->getVertex()->getLocation().getPt();
  const auto edge_vector = Normal(pt_1 - pt_0);
  const auto& face_plane = a_end->getFace()->getPlane();
  const auto& face_normal = face_plane.normal();
  // Normal normal_0 = getParaboloidSurfaceNormal(a_paraboloid, pt_0);
  // normal_0.approximatelyNormalize();
  // Normal normal_1 = getParaboloidSurfaceNormal(a_paraboloid, pt_1);
  // normal_1.approximatelyNormalize();
  Normal tgt_0 =
      computeTangentVectorAtPoint<ScalarType>(a_paraboloid, face_normal, pt_0);
  Normal tgt_1 =
      computeTangentVectorAtPoint<ScalarType>(a_paraboloid, face_normal, pt_1);
  const bool elliptic_face = a_paraboloid.a() * a_paraboloid.b() > ZERO &&
                             fabs(face_normal[2]) > MACHINE_EPSILON;

  if ((tgt_0[0] == ZERO && tgt_0[1] == ZERO && tgt_0[2] == ZERO) ||
      (tgt_1[0] == ZERO && tgt_1[1] == ZERO && tgt_1[2] == ZERO)) {
    *a_requires_nudge = true;
    return ReturnType::fromScalarConstant(ZERO);
  }

  if (!elliptic_face)  // The arc is from a hyperbola or parabola
  {
    const Normal n_cross_t0 = crossProduct(face_normal, tgt_0);
    // if (squaredMagnitude(edge_vector) < DISTANCE_EPSILON *
    // DISTANCE_EPSILON)
    // {
    //   if constexpr (!std::is_same<SurfaceOutputType,
    //   NoSurfaceOutput>::value)
    //   {
    //     auto surface_arc = RationalBezierArc(
    //         pt_1.toDoublePt(), 0.5 * (pt_0.toDoublePt() +
    //         pt_1.toDoublePt()), pt_0.toDoublePt(), 0.0);
    //     surface_arc.reset_start_point_id(
    //         reinterpret_cast<std::uintptr_t>(&pt_1));
    //     surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_0));
    //     a_surface->addArc(surface_arc);
    //   }
    // if (fabs(n_cross_t0 * tgt_1) < ANGLE_EPSILON &&
    //     fabs(ONE - fabs(dotProduct(edge_vector, tgt_0))) <
    //     ANGLE_EPSILON)
    //     {
    //   return ReturnType::fromScalarConstant(ZERO);
    // } else {
    //   *a_requires_nudge = true;
    //   return ReturnType::fromScalarConstant(ZERO);
    //   // }
    // } else {
    const ScalarType triple_prod = n_cross_t0 * tgt_1;
    if (fabs(triple_prod) < ANGLE_EPSILON) {  // The tangents are (close
                                              // to being) parallel
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        auto surface_arc = RationalBezierArc(
            pt_1.toDoublePt(), 0.5 * (pt_0.toDoublePt() + pt_1.toDoublePt()),
            pt_0.toDoublePt(), 0.0);
        surface_arc.reset_start_point_id(
            reinterpret_cast<std::uintptr_t>(&pt_1));
        surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_0));
        a_surface->addArc(surface_arc);
      }
      return ReturnType::fromScalarConstant(ZERO);
    } else {  // The tangents are NOT parallel
      // Compute control point
      const ScalarType lambda_1 = -(n_cross_t0 * edge_vector) / triple_prod;
      auto control_pt = Pt(pt_1 + lambda_1 * tgt_1);
      const auto pt0_to_cp = Normal(control_pt - pt_0);
      const ScalarType correction_normal = pt0_to_cp * face_normal;
      control_pt = control_pt - correction_normal * face_normal;

      // Orient tangents to point towards control point
      tgt_0 = (tgt_0 * pt0_to_cp) < ZERO ? -tgt_0 : tgt_0;
      tgt_1 = lambda_1 < ZERO ? -tgt_1 : tgt_1;

      // Construct Bezier arc from tangents
      const auto arc = RationalBezierArcBase<ScalarType>(
          pt_1, control_pt, pt_0, face_normal, a_paraboloid);
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        auto surface_arc = RationalBezierArc(
            pt_1.toDoublePt(), tgt_1.toDoubleNormal(), pt_0.toDoublePt(),
            tgt_0.toDoubleNormal(), face_normal.toDoubleNormal(),
            AlignedParaboloidBase<double>(a_paraboloid));
        surface_arc.reset_start_point_id(
            reinterpret_cast<std::uintptr_t>(&pt_1));
        surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_0));
        a_surface->addArc(surface_arc);
      }
      // if (arc.weight() < ZERO) {
      //   *a_requires_nudge = true;
      //   return ReturnType::fromScalarConstant(ZERO);
      // }
      return computeType3Contribution<ReturnType, ScalarType>(a_paraboloid,
                                                              arc);
    }
    // }
  } else  // The arc is from an ellipse
  {
    tgt_0.normalize();
    tgt_1.normalize();

    // Compute edge vectors and check if tangent is parallel to
    // edge
    Normal edge_0 = a_start->getVertex()->getLocation().getPt() -
                    a_start->getPreviousVertex()->getLocation().getPt();
    // if (squaredMagnitude(edge_0) < DISTANCE_EPSILON *
    // DISTANCE_EPSILON) {
    //   Normal edge_01 =
    //       a_start->getNextHalfEdge()->getVertex()->getLocation().getPt()
    //       - a_start->getVertex()->getLocation().getPt();
    //   edge_0 = (squaredMagnitude(edge_0) >
    //   squaredMagnitude(edge_01)) ? edge_0
    //                                                                   :
    //                                                                   edge_01;
    //   if (squaredMagnitude(edge_0) < DISTANCE_EPSILON *
    //   DISTANCE_EPSILON) {
    //     *a_requires_nudge = true;
    //     return ReturnType::fromScalarConstant(ZERO);
    //   }
    // }
    // edge_0.approximatelyNormalize();
    Normal edge_1 = a_end->getVertex()->getLocation().getPt() -
                    a_end->getPreviousVertex()->getLocation().getPt();
    // if (squaredMagnitude(edge_1) < DISTANCE_EPSILON *
    // DISTANCE_EPSILON) {
    //   Normal edge_11 =
    //       a_end->getNextHalfEdge()->getVertex()->getLocation().getPt()
    //       - a_end->getVertex()->getLocation().getPt();
    //   edge_1 = (squaredMagnitude(edge_1) >
    //   squaredMagnitude(edge_11)) ? edge_1
    //                                                                   :
    //                                                                   edge_11;
    //   if (squaredMagnitude(edge_1) < DISTANCE_EPSILON *
    //   DISTANCE_EPSILON) {
    //     *a_requires_nudge = true;
    //     return ReturnType::fromScalarConstant(ZERO);
    //   }
    // }
    // edge_1.approximatelyNormalize();

    // bool tgt_0_parallel_edge_0 = fabs(normal_0 * edge_0) <
    // 0.01; bool tgt_1_parallel_edge_1 = fabs(normal_1 * edge_1)
    // < 0.01;

    // if (!tgt_0_parallel_edge_0 &&
    //     !tgt_1_parallel_edge_1)  // We orient the tangent with
    //     the edge
    //                              // normal
    // {
    //   const Normal edge_normal_0 = crossProduct(face_normal,
    //   edge_0); const Normal edge_normal_1 =
    //   crossProduct(face_normal, edge_1); tgt_0 = (edge_normal_0
    //   * tgt_0 < ZERO) ? -tgt_0 : tgt_0; tgt_1 = (edge_normal_1
    //   * tgt_1 < ZERO) ? -tgt_1 : tgt_1;
    // } else {
    edge_0.normalize();
    edge_1.normalize();
    bool tgt_0_parallel_edge_0 =
        fabs(ONE - fabs(tgt_0 * edge_0)) < ANGLE_EPSILON;
    bool tgt_1_parallel_edge_1 =
        fabs(ONE - fabs(tgt_1 * edge_1)) < ANGLE_EPSILON;
    if (!tgt_0_parallel_edge_0 &&
        !tgt_1_parallel_edge_1)  // We orient the tangent with the
                                 // edge normal
    {
      const Normal edge_normal_0 = crossProduct(face_normal, edge_0);
      const Normal edge_normal_1 = crossProduct(face_normal, edge_1);
      tgt_0 = (edge_normal_0 * tgt_0 < ZERO) ? -tgt_0 : tgt_0;
      tgt_1 = (edge_normal_1 * tgt_1 < ZERO) ? -tgt_1 : tgt_1;
    } else {
      if constexpr (std::is_same_v<ScalarType, double>) {
        *a_requires_nudge = true;
        return ReturnType::fromScalarConstant(ZERO);
      }
      const Pt conic_center = conicCenter<ScalarType>(face_plane, a_paraboloid);
      if (squaredMagnitude(pt_0 - conic_center) <
              DISTANCE_EPSILON * DISTANCE_EPSILON ||
          squaredMagnitude(pt_1 - conic_center) <
              DISTANCE_EPSILON * DISTANCE_EPSILON) {
        return ReturnType::fromScalarConstant(ZERO);
      }
      if (squaredMagnitude(pt_0) < DISTANCE_EPSILON * DISTANCE_EPSILON &&
          squaredMagnitude(pt_1) < DISTANCE_EPSILON * DISTANCE_EPSILON &&
          face_plane.distance() < DISTANCE_EPSILON) {
        return ReturnType::fromScalarConstant(ZERO);
      }
      // Compute ellipse center an orient tangents accordingly
      auto center_to_pt_0 = Normal(pt_0 - conic_center);
      auto center_to_pt_1 = Normal(pt_1 - conic_center);
      center_to_pt_0.normalize();
      center_to_pt_1.normalize();
      Normal dummy_tgt_0 = crossProduct(face_normal, center_to_pt_0);
      Normal dummy_tgt_1 = crossProduct(face_normal, center_to_pt_1);
      assert(fabs(tgt_0 * dummy_tgt_0) > ANGLE_EPSILON);
      assert(fabs(tgt_1 * dummy_tgt_1) > ANGLE_EPSILON);
      tgt_0 = (tgt_0 * dummy_tgt_0) < ZERO ? -tgt_0 : tgt_0;
      tgt_1 = (tgt_1 * dummy_tgt_1) > ZERO ? -tgt_1 : tgt_1;
      // At this point, the tangents form a valid arc (but they
      // may be oriented in the wrong direction)
      if (!tgt_0_parallel_edge_0) {
        const Normal edge_normal_0 = crossProduct(face_normal, edge_0);
        if (edge_normal_0 * tgt_0 < ZERO) {
          tgt_0 = -tgt_0;
          tgt_1 = -tgt_1;
        }
      } else if (!tgt_1_parallel_edge_1) {
        const Normal edge_normal_1 = crossProduct(face_normal, edge_1);
        if (edge_normal_1 * tgt_1 < ZERO) {
          tgt_0 = -tgt_0;
          tgt_1 = -tgt_1;
        }
      } else {
        if (a_paraboloid.a() < ZERO) {
          tgt_0 = -tgt_0;
          tgt_1 = -tgt_1;
        }
      }
    }
    // }
    UnsignedIndex_t split_counter = 0;
    return computeType3ContributionWithSplit<ReturnType, ScalarType>(
        a_paraboloid, face_normal, pt_1, pt_1, pt_0, tgt_1, tgt_0,
        a_requires_nudge, &split_counter, a_surface);
  }

  return ReturnType::fromScalarConstant(ZERO);
}  // namespace IRL

// TODO
template <class ReturnType, class ScalarType, class HalfEdgeType,
          class SurfaceOutputType>
ReturnType orientAndApplyType3CorrectionWithGradients(
    const AlignedParaboloid& a_paraboloid, HalfEdgeType* a_start,
    HalfEdgeType* a_end, SurfaceOutputType* a_surface, bool* a_requires_nudge) {
  /* Defining constants and types */
  using gradient_type = typename ReturnType::gradient_type;
  using pt_type = PtWithGradient<gradient_type>;
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  using Plane = PlaneBase<ScalarType>;
  using AlignedParaboloid = AlignedParaboloidBase<ScalarType>;
  const ScalarType EPSILON =
      maximum(machine_epsilon<ScalarType>(), static_cast<ScalarType>(1.0e-24q));
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType FOUR = static_cast<ScalarType>(4);
  const ScalarType HALF = ONE / TWO;
  const ScalarType ONEQUARTER = ONE / FOUR;
  const ScalarType THREEQUARTERS = HALF + ONEQUARTER;
  const ScalarType TEN = static_cast<ScalarType>(10);
  const ScalarType ONE_HUNDRED = static_cast<ScalarType>(100);

  /* Function */
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
  face_normal_grad[0].setGradRx(ZERO);
  face_normal_grad[1].setGradRx(face_normal[2]);
  face_normal_grad[2].setGradRx(-face_normal[1]);
  face_normal_grad[0].setGradRy(-face_normal[2]);
  face_normal_grad[1].setGradRy(ZERO);
  face_normal_grad[2].setGradRy(face_normal[0]);
  face_normal_grad[0].setGradRz(face_normal[1]);
  face_normal_grad[1].setGradRz(-face_normal[0]);
  face_normal_grad[2].setGradRz(ZERO);

  // Compute tangents and their gradient
  auto tgt_withgrad_0 =
      computeTangentVectorAndGradientAtPoint<ScalarType, pt_type>(
          a_paraboloid, face_normal_withgrad, pt_withgrad_0);
  auto tgt_withgrad_1 =
      computeTangentVectorAndGradientAtPoint<ScalarType, pt_type>(
          a_paraboloid, face_normal_withgrad, pt_withgrad_1);
  auto tgt_0 = Normal::fromPt(tgt_withgrad_0.getPt());
  auto tgt_1 = Normal::fromPt(tgt_withgrad_1.getPt());
  const bool elliptic_face =
      a_paraboloid.a() * a_paraboloid.b() > ZERO && face_normal[2] != ZERO;

  if (!elliptic_face)  // The arc is from a hyperbola or parabola
  {
    const Normal n_cross_t0 = crossProduct(face_normal, tgt_0);
    if (fabs(n_cross_t0 * tgt_1) <
        ONE_HUNDRED * EPSILON)  // The tangents are (close to
                                // being) parallel
    {
      // We consider an arc with zero contribution
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        auto surface_arc = RationalBezierArc(
            pt_1.toDoublePt(), 0.5 * (pt_0.toDoublePt() + pt_1.toDoublePt()),
            pt_0.toDoublePt(), 0.0);
        surface_arc.reset_start_point_id(
            reinterpret_cast<std::uintptr_t>(&pt_1));
        surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_0));
        a_surface->addArc(surface_arc);
      }
      return ReturnType::fromScalarConstant(ZERO);
    } else {  // The tangents are NOT parallel
      // Compute control point
      const ScalarType lambda_1 =
          -(n_cross_t0 * edge_vector) / safelyEpsilon(n_cross_t0 * tgt_1);
      const auto control_pt = Pt(pt_1 + lambda_1 * tgt_1);
      // Orient tangents to point towards control point
      tgt_withgrad_0 = (tgt_0 * Normal(control_pt - pt_0)) < ZERO
                           ? -tgt_withgrad_0
                           : tgt_withgrad_0;
      tgt_withgrad_1 = (tgt_1 * Normal(control_pt - pt_1)) < ZERO
                           ? -tgt_withgrad_1
                           : tgt_withgrad_1;
      // Construct Bezier arc from tangents
      const auto arc_with_gradient =
          RationalBezierArcWithGradientBase<pt_type, ScalarType>(
              pt_withgrad_1, tgt_withgrad_1, pt_withgrad_0, tgt_withgrad_0,
              face_normal_withgrad, a_paraboloid);
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        auto surface_arc = RationalBezierArc(
            pt_1.toDoublePt(),
            arc_with_gradient.control_point().getPt().toDoublePt(),
            pt_0.toDoublePt(), static_cast<double>(arc_with_gradient.weight()));
        surface_arc.reset_start_point_id(
            reinterpret_cast<std::uintptr_t>(&pt_1));
        surface_arc.reset_end_point_id(reinterpret_cast<std::uintptr_t>(&pt_0));
        a_surface->addArc(surface_arc);
      }
      return computeType3ContributionWithGradient<ReturnType, ScalarType>(
          a_paraboloid, arc_with_gradient);
    }
  } else  // The arc is from an ellipse
  {
    // Compute edge vectors and check if tangent is parallel to
    // edge
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
        fabs(ONE - fabs(tgt_0 * edge_0)) < ONE_HUNDRED * EPSILON;
    bool tgt_1_parallel_edge_1 =
        fabs(ONE - fabs(tgt_1 * edge_1)) < ONE_HUNDRED * EPSILON;
    if (!tgt_0_parallel_edge_0 &&
        !tgt_1_parallel_edge_1)  // We orient the tangent with the
                                 // edge normal
    {
      const Normal edge_normal_0 = crossProduct(face_normal, edge_0);
      const Normal edge_normal_1 = crossProduct(face_normal, edge_1);
      tgt_withgrad_0 =
          (edge_normal_0 * tgt_0 < ZERO) ? -tgt_withgrad_0 : tgt_withgrad_0;
      tgt_withgrad_1 =
          (edge_normal_1 * tgt_1 < ZERO) ? -tgt_withgrad_1 : tgt_withgrad_1;
    } else {
      // Compute ellipse center an orient tangents accordingly
      const Pt conic_center = conicCenter(face_plane, a_paraboloid);
      auto center_to_pt_0 = Normal(pt_0 - conic_center);
      center_to_pt_0.normalize();
      auto center_to_pt_1 = Normal(pt_1 - conic_center);
      center_to_pt_1.normalize();
      Normal dummy_tgt_0 = crossProduct(face_normal, center_to_pt_0);
      Normal dummy_tgt_1 = crossProduct(face_normal, center_to_pt_1);
      assert(fabs(tgt_0 * dummy_tgt_0) > TEN * EPSILON);
      assert(fabs(tgt_1 * dummy_tgt_1) > TEN * EPSILON);
      tgt_0 = (tgt_0 * dummy_tgt_0) < ZERO ? -tgt_0 : tgt_0;
      tgt_withgrad_0 =
          (tgt_0 * dummy_tgt_0) < ZERO ? -tgt_withgrad_0 : tgt_withgrad_0;
      tgt_1 = (tgt_1 * dummy_tgt_1) > ZERO ? -tgt_1 : tgt_1;
      tgt_withgrad_1 =
          (tgt_1 * dummy_tgt_1) > ZERO ? -tgt_withgrad_1 : tgt_withgrad_1;
      // At this point, the tangents form a valid arc (but they
      // may be oriented in the wrong direction)
      if (!tgt_0_parallel_edge_0) {
        const Normal edge_normal_0 = crossProduct(face_normal, edge_0);
        if (edge_normal_0 * tgt_0 < ZERO) {
          tgt_withgrad_0 = -tgt_withgrad_0;
          tgt_withgrad_1 = -tgt_withgrad_1;
        }
      } else if (!tgt_1_parallel_edge_1) {
        const Normal edge_normal_1 = crossProduct(face_normal, edge_1);
        if (edge_normal_1 * tgt_1 < ZERO) {
          tgt_withgrad_0 = -tgt_withgrad_0;
          tgt_withgrad_1 = -tgt_withgrad_1;
        }
      } else {
        // The arcs should have been sorted in this case!! This
        // allows to project the mid-point onto the arc, and check
        // that it belongs to the face
        const Pt avg_pt = HALF * (pt_0 + pt_1);
        auto avg_tgt = Normal(HALF * (tgt_0 + tgt_1));
        if (squaredMagnitude(avg_tgt) <
            static_cast<ScalarType>(1.0e6) * EPSILON * EPSILON) {
          avg_tgt = Normal(ONEQUARTER * tgt_0 + THREEQUARTERS * tgt_1);
        }
        avg_tgt.approximatelyNormalize();
        Pt proj_test = projectPtAlongHalfLineOntoParaboloid<ScalarType>(
            a_paraboloid, avg_tgt, avg_pt);
        // Check if test projected point is inside the face
        UnsignedIndex_t best_proj_dir = 0;
        if (fabs(face_normal[best_proj_dir]) < fabs(face_normal[1]))
          best_proj_dir = 1;
        if (fabs(face_normal[best_proj_dir]) < fabs(face_normal[2]))
          best_proj_dir = 2;
        const UnsignedIndex_t start_id_for_search = (best_proj_dir + 1) % 3;
        auto current_half_edge = a_start;
        bool pt_internal_to_polygon = false;
        do {
          const Pt& location_0 =
              current_half_edge->getPreviousVertex()->getLocation().getPt();
          const Pt& location_1 =
              current_half_edge->getVertex()->getLocation().getPt();
          if (isPtBeforeIntersectionWithEdgeWithComponent<ScalarType>(
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
    return computeType3ContributionWithGradientWithSplit<ReturnType,
                                                         ScalarType>(
        a_paraboloid, face_normal_withgrad, pt_withgrad_1, pt_withgrad_1,
        pt_withgrad_0, tgt_withgrad_1, tgt_withgrad_0, a_surface);
  }

  return ReturnType::fromScalarConstant(ZERO);
}  // namespace IRL

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
      // New vertex from intersection, remove it and patch
      // half-edges
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

template <class ScalarType>
AlignedParaboloidBase<Quad_t> nudgeParaboloid(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const UnsignedIndex_t a_nudge_iter) {
  std::random_device rd;
  std::mt19937 gen(a_nudge_iter);
  std::uniform_real_distribution distr(-1.0, 1.0);

  const Quad_t nudge_epsilon = 0.0q;  // 1.0e6q * distance_epsilon<Quad_t>();
  return AlignedParaboloidBase<Quad_t>(
      {static_cast<Quad_t>(a_paraboloid.a() + nudge_epsilon * distr(gen)),
       static_cast<Quad_t>(a_paraboloid.b() + nudge_epsilon * distr(gen))});
}

/********************** Move poly to find stable topology ****************/
template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType,
          class SurfaceOutputType>
void nudgePolyhedron(SegmentedHalfEdgePolyhedronType* a_polytope,
                     HalfEdgePolytopeType* a_complete_polytope,
                     const UnsignedIndex_t a_nudge_iter,
                     SurfaceOutputType* a_surface) {
  using pt_type =
      typename SegmentedHalfEdgePolyhedronType::vertex_type::pt_type;
  static_assert(std::is_same_v<typename pt_type::value_type, Quad_t>);

  std::random_device rd;
  std::mt19937 gen(a_nudge_iter);
  std::uniform_real_distribution distr(-1.0, 1.0);

  const Quad_t nudge_epsilon = 1.0e8q * distance_epsilon<Quad_t>();

  // std::cout << " NUDGE " << a_nudge_iter << std::endl;

  auto center = a_polytope->calculateCentroid();
  auto converted_center = PtBase<Quad_t>(static_cast<Quad_t>(center[0]),
                                         static_cast<Quad_t>(center[1]),
                                         static_cast<Quad_t>(center[2]));
  const auto nudge_translation =
      nudge_epsilon * NormalBase<Quad_t>(static_cast<Quad_t>(distr(gen)),
                                         static_cast<Quad_t>(distr(gen)),
                                         static_cast<Quad_t>(distr(gen)));
  const auto nudge_rotation =
      2.0q * nudge_epsilon *
      NormalBase<Quad_t>(static_cast<Quad_t>(distr(gen)),
                         static_cast<Quad_t>(distr(gen)),
                         static_cast<Quad_t>(distr(gen)));

  ReferenceFrameBase<Quad_t> frame(NormalBase<Quad_t>(1, 0, 0),
                                   NormalBase<Quad_t>(0, 1, 0),
                                   NormalBase<Quad_t>(0, 0, 1));
  UnitQuaternionBase<Quad_t> x_rotation(nudge_rotation[0], frame[0]);
  UnitQuaternionBase<Quad_t> y_rotation(nudge_rotation[1], frame[1]);
  UnitQuaternionBase<Quad_t> z_rotation(nudge_rotation[2], frame[2]);
  auto total_rotation = x_rotation * y_rotation * z_rotation;
  total_rotation.normalize();
  frame = total_rotation * frame;

  for (UnsignedIndex_t v = 0; v < a_polytope->getNumberOfVertices(); ++v) {
    const PtBase<Quad_t> original_pt =
        a_polytope->getVertex(v)->getLocation().getPt() - converted_center;
    PtBase<Quad_t> projected_location(0, 0, 0);
    auto& pt = projected_location.getPt();
    for (UnsignedIndex_t n = 0; n < 3; ++n) {
      pt[n] = frame[n] * original_pt;
    }
    pt += converted_center + nudge_translation;
    a_polytope->getVertex(v)->setLocation(projected_location);
  }

  if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
    a_surface->clear();
  }
}

template <class DoubleSegmentedHalfEdgePolytopeType,
          class DoubleHalfEdgePolytopeType, class QuadHalfEdgePolytopeType>
void convertPolytopeFromDoubleToQuadPrecision(
    DoubleSegmentedHalfEdgePolytopeType* a_polytope,
    DoubleHalfEdgePolytopeType* a_complete_polytope,
    QuadHalfEdgePolytopeType* a_converted_polytope) {
  /* Convert polytope type */
  using converted_pt_type = PtBase<Quad_t>;
  using converted_vertex_type = VertexParaboloid<converted_pt_type>;
  using converted_halfedge_type = HalfEdgeParaboloid<converted_vertex_type>;
  using converted_face_type = FaceParaboloid<converted_halfedge_type>;
  const UnsignedIndex_t converted_kMaxHalfEdges =
      DoubleHalfEdgePolytopeType::maxHalfEdges;
  const UnsignedIndex_t converted_kMaxVertices =
      DoubleHalfEdgePolytopeType::maxVertices;
  const UnsignedIndex_t converted_kMaxFaces =
      DoubleHalfEdgePolytopeType::maxFaces;

  // /* First, triangulate faces */
  // triangulatePolytope(a_polytope, a_complete_polytope);

  /* Calculate number of half edges, vertices and faces*/
  UnsignedIndex_t number_of_half_edges = 0;
  for (const auto& face : (*a_polytope)) {
    const auto starting_half_edge = face->getStartingHalfEdge();
    auto current_half_edge = starting_half_edge;
    do {
      number_of_half_edges++;
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != starting_half_edge);
  }
  const UnsignedIndex_t number_of_vertices = a_polytope->getNumberOfVertices();
  const UnsignedIndex_t number_of_faces = a_polytope->getNumberOfFaces();

  /* Convert vertices */
  std::vector<PtBase<Quad_t>> pt_list;
  pt_list.resize(number_of_vertices);
  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    const auto old_pt = a_polytope->getVertex(v)->getLocation();
    pt_list[v] = PtBase<Quad_t>(static_cast<Quad_t>(old_pt[0]),
                                static_cast<Quad_t>(old_pt[1]),
                                static_cast<Quad_t>(old_pt[2]));
  }

  /* Create face -> halfedge and halfedge -> opposite halfedge
   * mapping */
  std::vector<std::vector<UnsignedIndex_t>> face_mapping;
  face_mapping.resize(number_of_faces);

  UnsignedIndex_t f = 0;
  for (const auto& face : (*a_polytope)) {
    const auto starting_half_edge = face->getStartingHalfEdge();
    auto current_half_edge = starting_half_edge;
    do {
      bool found = false;
      for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
        if (current_half_edge->getVertex() == a_polytope->getVertex(v)) {
          face_mapping[f].push_back(v);
          found = true;
          break;
        }
      }
      assert(found);
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != starting_half_edge);
    f++;
  }

  /* Create converted polytope */
  BREPToHalfEdge<converted_pt_type, converted_vertex_type,
                 converted_halfedge_type, converted_face_type,
                 converted_kMaxHalfEdges, converted_kMaxVertices,
                 converted_kMaxFaces>::setHalfEdgeVersion(pt_list, face_mapping,
                                                          a_converted_polytope);
}

template <class HalfEdgeType>
HalfEdgeType* getNextNonIntersection(HalfEdgeType* a_half_edge) {
  if (a_half_edge->getVertex()->needsToSeek()) {
    return a_half_edge;
  } else if (a_half_edge->getNextHalfEdge()->getVertex()->needsToSeek()) {
    return a_half_edge->getNextHalfEdge();
  } else {
    return a_half_edge->getNextHalfEdge()->getNextHalfEdge();
  }
}

template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType,
          class AlignedParaboloidType, class ScalarType>
void triangulatePolytope(SegmentedHalfEdgePolyhedronType* a_polytope,
                         HalfEdgePolytopeType* a_complete_polytope,
                         const AlignedParaboloidType& a_aligned_paraboloid,
                         const ScalarType a_nudge_epsilon,
                         bool* a_requires_nudge) {
  using VertexType = typename HalfEdgePolytopeType::vertex_type;
  using PtType = typename VertexType::pt_type;
  using HalfEdgeType = typename HalfEdgePolytopeType::half_edge_type;
  using FaceType = typename HalfEdgePolytopeType::face_type;

  // // Loop over all faces. Completely face independent procedure
  // const auto nfaces = a_polytope->getNumberOfFaces();
  // for (UnsignedIndex_t f = 0; f < nfaces; ++f) {
  //   auto face = (*a_polytope)[f];
  //   auto starting_half_edge = face->getStartingHalfEdge();
  //   auto nf_start = starting_half_edge->getNextHalfEdge();  //
  //   New face start auto next = nf_start->getNextHalfEdge();
  //   auto next_base_edge = next->getNextHalfEdge(); while
  //   (next->getNextHalfEdge() != starting_half_edge) {
  //     // Create new triangular face
  //     auto new_face =
  //     a_complete_polytope->getNewFace(FaceType(nf_start));
  //     a_polytope->addFace(new_face);

  //     // Creating new half-edge to close new triangular face
  //     auto new_half_edge =
  //     a_complete_polytope->getNewHalfEdge(HalfEdgeType(
  //         starting_half_edge->getVertex(), next, nf_start,
  //         new_face));
  //     nf_start->setPreviousHalfEdge(new_half_edge);
  //     nf_start->setFace(new_face);
  //     next->setNextHalfEdge(new_half_edge);
  //     next->setFace(new_face);

  //     // Creating new opposite half-edge
  //     // Use of starting_half_edge was wrong here if not a face
  //     with 4
  //     // vertices
  //     auto new_opposite_half_edge =
  //     a_complete_polytope->getNewHalfEdge(
  //         HalfEdgeType(next->getVertex(), nullptr,
  //         next_base_edge, nullptr));
  //     setMutualOpposites(new_half_edge,
  //     new_opposite_half_edge);

  //     // Update for next iteration
  //     nf_start = new_opposite_half_edge;
  //     next = next_base_edge;
  //     next_base_edge = next->getNextHalfEdge();
  //     assert(starting_half_edge->getVertex()->checkValidHalfEdgeCycle());
  //     assert(new_half_edge->getVertex()->checkValidHalfEdgeCycle());
  //   }
  //   nf_start->setPreviousHalfEdge(starting_half_edge);
  //   starting_half_edge->setNextHalfEdge(nf_start);
  //   nf_start->setFace(face);
  // }
  // assert(a_polytope->checkValidHalfEdgeStructure());

  //   // Loop over all faces. Completely face independent procedure
  // const auto nfaces = a_polytope->getNumberOfFaces();
  // UnsignedIndex_t new_faces = 0;
  // for (UnsignedIndex_t f = 0; f < nfaces; ++f) {
  //   auto face = (*a_polytope)[f];
  //   auto starting_half_edge = face->getStartingHalfEdge();
  //   auto half_edge = starting_half_edge->getNextHalfEdge();
  //   auto next = half_edge->getNextHalfEdge();
  //   if (half_edge == starting_half_edge || next == starting_half_edge) {
  //     continue;
  //   }
  //   auto next_next = next->getNextHalfEdge();
  //   while (next_next != starting_half_edge) {
  //     // Create new triangular face
  //     a_polytope->addFace(a_complete_polytope->getNewFace(FaceType(half_edge)));
  //     auto new_face = (*a_polytope)[nfaces + new_faces++];
  //     // Creating new half-edge to close new triangular face
  //     auto new_half_edge =
  //     a_complete_polytope->getNewHalfEdge(HalfEdgeType(
  //         starting_half_edge->getVertex(), next, half_edge, new_face));
  //     half_edge->setPreviousHalfEdge(new_half_edge);
  //     half_edge->setFace(new_face);
  //     next->setNextHalfEdge(new_half_edge);
  //     next->setFace(new_face);

  //     // Creating new opposite half-edge
  //     auto new_opposite_half_edge = a_complete_polytope->getNewHalfEdge(
  //         HalfEdgeType(next->getVertex(), starting_half_edge, next_next,
  //         face));
  //     setMutualOpposites(new_half_edge, new_opposite_half_edge);
  //     starting_half_edge->setNextHalfEdge(new_opposite_half_edge);
  //     next_next->setPreviousHalfEdge(new_opposite_half_edge);

  //     // Closing old face
  //     half_edge = new_opposite_half_edge;
  //     next = next_next;
  //     next_next = next->getNextHalfEdge();

  //     assert(starting_half_edge->getVertex()->checkValidHalfEdgeCycle());
  //     assert(new_half_edge->getVertex()->checkValidHalfEdgeCycle());
  //   }
  // }

  // Loop over all faces. Completely face independent procedure
  const auto nfaces = a_polytope->getNumberOfFaces();
  UnsignedIndex_t new_faces = 0;
  for (UnsignedIndex_t f = 0; f < nfaces; ++f) {
    auto face = (*a_polytope)[f];
    /* Only triangulate faces with intersections */
    if (face->getNumberOfIntersections() > 0) {
      auto starting_half_edge =
          getNextNonIntersection(face->getStartingHalfEdge());
      auto half_edge0 = starting_half_edge->getNextHalfEdge();
      auto half_edge1 = getNextNonIntersection(half_edge0);
      auto next0 = half_edge1->getNextHalfEdge();
      auto next1 = getNextNonIntersection(next0);
      if (half_edge1 == starting_half_edge || next1 == starting_half_edge) {
        continue;
      }
      auto next_next0 = next1->getNextHalfEdge();
      auto next_next1 = getNextNonIntersection(next_next0);

      while (next_next1 != starting_half_edge) {
        // Create new triangular face
        a_polytope->addFace(
            a_complete_polytope->getNewFace(FaceType(half_edge1)));
        auto new_face = (*a_polytope)[nfaces + new_faces++];
        // Creating new half-edge to close new triangular face
        auto new_half_edge = a_complete_polytope->getNewHalfEdge(HalfEdgeType(
            starting_half_edge->getVertex(), next1, half_edge0, new_face));
        new_face->setStartingHalfEdge(new_half_edge);
        half_edge0->setPreviousHalfEdge(new_half_edge);
        half_edge0->setFace(new_face);
        if (half_edge0 != half_edge1) {
          half_edge1->setFace(new_face);
          if (half_edge0->getNextHalfEdge() != half_edge1) {
            half_edge0->getNextHalfEdge()->setFace(new_face);
          }
        }
        next1->setNextHalfEdge(new_half_edge);
        next0->setFace(new_face);
        if (next0 != next1) {
          next1->setFace(new_face);
          if (next0->getNextHalfEdge() != next1) {
            next0->getNextHalfEdge()->setFace(new_face);
          }
        }

        // Creating new opposite half-edge
        auto new_opposite_half_edge =
            a_complete_polytope->getNewHalfEdge(HalfEdgeType(
                next1->getVertex(), starting_half_edge, next_next0, face));
        setMutualOpposites(new_half_edge, new_opposite_half_edge);
        starting_half_edge->setNextHalfEdge(new_opposite_half_edge);
        next_next0->setPreviousHalfEdge(new_opposite_half_edge);

        // Closing old face
        half_edge0 = new_opposite_half_edge;
        half_edge1 = half_edge0;
        next0 = next_next0;
        next1 = next_next1;
        next_next0 = next1->getNextHalfEdge();
        next_next1 = getNextNonIntersection(next_next0);

        assert(starting_half_edge->getVertex()->checkValidHalfEdgeCycle());
        assert(new_half_edge->getVertex()->checkValidHalfEdgeCycle());
      }
    }
  }

  if (new_faces > 0) {
    /* Compute intersections on new edges */
    const bool elliptic = a_aligned_paraboloid.a() * a_aligned_paraboloid.b() >
                          static_cast<ScalarType>(0);
    const ScalarType nudge_epsilon_sq = a_nudge_epsilon * a_nudge_epsilon;
    StackVector<PtType, 2> edge_intercepts;
    for (UnsignedIndex_t f = nfaces; f < nfaces + new_faces; ++f) {
      auto face = (*a_polytope)[f];
      face->markAsVisited();
      auto current_edge = face->getStartingHalfEdge();
      const auto& vertex_start = current_edge->getPreviousVertex();
      const auto& vertex_end = current_edge->getVertex();
      if (elliptic) {
        if (a_aligned_paraboloid.a() > static_cast<ScalarType>(0)) {
          if (vertex_start->isNotClipped() && vertex_end->isNotClipped()) {
            continue;
          }
        } else if (vertex_start->isClipped() && vertex_end->isClipped()) {
          continue;
        }
      }
      const auto& edge_start = vertex_start->getLocation();
      const auto& edge_end = vertex_end->getLocation();
      checkAndFindIntercepts<PtType, ScalarType>(
          a_aligned_paraboloid, edge_start, edge_end, &edge_intercepts,
          a_nudge_epsilon, elliptic);
      if (edge_intercepts.size() == 1) {
        if (squaredMagnitude(edge_intercepts[0] - edge_start) <
                nudge_epsilon_sq ||
            squaredMagnitude(edge_intercepts[0] - edge_end) <
                nudge_epsilon_sq) {
          *a_requires_nudge = true;
          break;
        }
        placeSingleIntercept(edge_intercepts[0], current_edge, a_polytope,
                             a_complete_polytope);
        current_edge->getPreviousVertex()->setAsUnnecessaryToSeek();
      } else if (edge_intercepts.size() == 2) {
        if (squaredMagnitude(edge_intercepts[0] - edge_start) <
                nudge_epsilon_sq ||
            squaredMagnitude(edge_intercepts[1] - edge_end) <
                nudge_epsilon_sq ||
            squaredMagnitude(edge_intercepts[0] - edge_intercepts[1]) <
                nudge_epsilon_sq) {
          *a_requires_nudge = true;
          break;
        }
        placeDoubleIntercept(edge_intercepts, current_edge, a_polytope,
                             a_complete_polytope);
        current_edge->getPreviousVertex()->setAsUnnecessaryToSeek();
        current_edge->getPreviousHalfEdge()
            ->getPreviousVertex()
            ->setAsUnnecessaryToSeek();
      }
    }

    // if (!a_polytope->checkValidHalfEdgeStructure()) {
    //   std::cout << (*a_polytope) << std::endl;
    //   // std::cout << (*a_complete_polytope) << std::endl;
    // }

    assert(a_polytope->checkValidHalfEdgeStructure());

    if (!(*a_requires_nudge)) {
      /* Make starting H-edge of face an entry and update number of
       * intersections
       */
      for (UnsignedIndex_t f = 0; f < nfaces + new_faces; ++f) {
        auto face = (*a_polytope)[f];
        if (face->getNumberOfIntersections() > 0 || f >= nfaces) {
          face->clearIntersections();
          const auto starting_half_edge = face->getStartingHalfEdge();
          auto current_edge = starting_half_edge;
          do {
            if (current_edge->getVertex()->doesNotNeedToSeek()) {
              face->addIntersection();
            }
            current_edge = current_edge->getNextHalfEdge();
          } while (current_edge != starting_half_edge);
        }
      }
      for (UnsignedIndex_t f = nfaces; f < nfaces + new_faces; ++f) {
        auto face = (*a_polytope)[f];
        const auto starting_half_edge = face->getStartingHalfEdge();
        auto current_edge = starting_half_edge;
        do {
          if (current_edge->getVertex()->doesNotNeedToSeek()) {
            if (current_edge->getPreviousVertex()->isClipped()) {
              face->setStartingHalfEdge(current_edge);
              break;
            }
          } else if (current_edge->getVertex()->isNotClipped() &&
                     current_edge->getPreviousVertex()->doesNotNeedToSeek()) {
            face->setStartingHalfEdge(current_edge->getPreviousHalfEdge());
            break;
          }
          current_edge = current_edge->getNextHalfEdge();
        } while (current_edge != starting_half_edge);
      }
    }
  }
}

// Assumes paraboloid of function 0 = a*x^2 + b*y^2 + z.
// We will truncate the polyhedron to exist in the region
// q < a*x^2 + b*y^2 + z
template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class AligneParaboloidType,
          class SurfaceOutputType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
formParaboloidIntersectionBases(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AligneParaboloidType& a_aligned_paraboloid,
    const UnsignedIndex_t a_nudge_iter, SurfaceOutputType* a_surface) {
  using vertex_type = typename SegmentedHalfEdgePolyhedronType::vertex_type;
  using pt_type = typename vertex_type::pt_type;
  using face_type = typename SegmentedHalfEdgePolyhedronType::face_type;
  using half_edge_type = typename HalfEdgePolytopeType::half_edge_type;

  /* Defining type aliases (needed to ensure precision is
   * consistent)*/
  using ScalarType = typename pt_type::value_type;
  static_assert(std::is_same_v<typename pt_type::value_type, ScalarType>);
  static_assert(
      std::is_same_v<typename half_edge_type::value_type, ScalarType>);
  static_assert(std::is_same_v<typename face_type::value_type, ScalarType>);
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  using Plane = PlaneBase<ScalarType>;
  using AlignedParaboloid = AlignedParaboloidBase<ScalarType>;

  using converted_pt_type = PtBase<Quad_t>;
  using converted_vertex_type = VertexParaboloid<converted_pt_type>;
  using converted_halfedge_type = HalfEdgeParaboloid<converted_vertex_type>;
  using converted_face_type = FaceParaboloid<converted_halfedge_type>;
  const UnsignedIndex_t converted_kMaxHalfEdges =
      HalfEdgePolytopeType::maxHalfEdges;
  const UnsignedIndex_t converted_kMaxVertices =
      HalfEdgePolytopeType::maxVertices;
  const UnsignedIndex_t converted_kMaxFaces = HalfEdgePolytopeType::maxFaces;
  using converted_complete_polytope_type =
      HalfEdgePolyhedronParaboloid<converted_pt_type, converted_vertex_type,
                                   converted_halfedge_type, converted_face_type,
                                   converted_kMaxHalfEdges,
                                   converted_kMaxVertices, converted_kMaxFaces>;

  /* Defining constants */
  const ScalarType MACHINE_EPSILON = machine_epsilon<ScalarType>();
  const ScalarType DISTANCE_EPSILON = distance_epsilon<ScalarType>();
  const ScalarType ANGLE_EPSILON = angle_epsilon<ScalarType>();
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType THREE = static_cast<ScalarType>(3);
  const ScalarType FIVE = static_cast<ScalarType>(5);
  const ScalarType TEN = static_cast<ScalarType>(10);
  const ScalarType ONE_HUNDRED = static_cast<ScalarType>(100);

  assert(!(a_surface != nullptr &&
           std::is_same<SurfaceOutputType, NoSurfaceOutput>::value));

  ReturnType full_moments = ReturnType::fromScalarConstant(ZERO);

  // Initialising variables for handling degenerate cases
  bool requires_nudge = false;
  const ScalarType nudge_epsilon = DISTANCE_EPSILON;
  const ScalarType nudge_epsilon_sq = nudge_epsilon * nudge_epsilon;

  if (a_nudge_iter >= 100) {
    std::cout << "Nudged more than 100 times " << std::endl;
    return ReturnType::fromScalarConstant(-static_cast<ScalarType>(9999999.9));
  }

  // /* First, triangulate faces */
  // if (a_nudge_iter == 0) {
  //   triangulatePolytope(a_polytope, a_complete_polytope);
  // }

  // // Recalculate in face plane information
  // for (auto& face : (*a_polytope)) {
  //   const auto starting_half_edge = face->getStartingHalfEdge();
  //   const auto& p0 =
  //       starting_half_edge->getPreviousVertex()->getLocation().getPt();
  //   const auto& p1 =
  //   starting_half_edge->getVertex()->getLocation().getPt(); const auto&
  //   p2 = starting_half_edge->getNextHalfEdge()
  //                        ->getVertex()
  //                        ->getLocation()
  //                        .getPt();
  //   auto normal = Normal(
  //       (p1[1] - p0[1]) * (p2[2] - p1[2]) - (p1[2] - p0[2]) * (p2[1] -
  //       p1[1]), (p1[2] - p0[2]) * (p2[0] - p1[0]) - (p1[0] - p0[0]) *
  //       (p2[2]
  //       - p1[2]), (p1[0] - p0[0]) * (p2[1] - p1[1]) - (p1[1] - p0[1]) *
  //       (p2[0] - p1[0]));
  //   ScalarType squaredMag = squaredMagnitude(normal);
  //   if (squaredMag < DISTANCE_EPSILON * DISTANCE_EPSILON) {
  //     normal = Normal(ZERO, ZERO, ZERO);
  //   } else {
  //     normal /= sqrt(squaredMag);
  //   }
  //   face->setPlane(Plane(normal, normal * p0));
  // }

  // Identify elliptic case
  const bool elliptic =
      a_aligned_paraboloid.a() * a_aligned_paraboloid.b() > ZERO;

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
  if (elliptic && a_aligned_paraboloid.a() > ZERO &&
      number_of_vertices_above == 0) {
    // This is needed to that we can reset the polytope after the
    // return
    for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
      auto& vertex = *(a_polytope->getVertex(v));
      vertex.setToSeek();
    }
    // Whole volume below
    return ReturnType::calculateMoments(a_polytope);
  }

  if (elliptic && a_aligned_paraboloid.a() < ZERO &&
      number_of_vertices_above == starting_number_of_vertices) {
    // This is needed to that we can reset the polytope after the
    // return
    for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
      auto& vertex = *(a_polytope->getVertex(v));
      vertex.setToSeek();
    }
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
  //      !elliptic || (elliptic && a_aligned_paraboloid.a() >
  //      0.0);
  const bool check_from_unclipped = true;
  //      !elliptic || (elliptic && !check_from_clipped);

  // Compute gradients of polyhedron corners (if gradients are
  // requested)
  if constexpr (has_embedded_gradient<ReturnType>::value) {
    for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
      auto& pt = a_polytope->getVertex(v)->getLocation().getPt();
      auto& gradient = a_polytope->getVertex(v)->getLocation().getData();
      gradient[0].setGradTx(-ONE);
      gradient[1].setGradTy(-ONE);
      gradient[2].setGradTz(-ONE);
      gradient[0].setGradRx(ZERO);
      gradient[1].setGradRx(pt[2]);
      gradient[2].setGradRx(-pt[1]);
      gradient[0].setGradRy(-pt[2]);
      gradient[1].setGradRy(ZERO);
      gradient[2].setGradRy(pt[0]);
      gradient[0].setGradRz(pt[1]);
      gradient[1].setGradRz(-pt[0]);
      gradient[2].setGradRz(ZERO);
    }
  }

  // Compute intersections
  StackVector<pt_type, 2> edge_intercepts;
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    vertex.setToSeek();
    if (check_from_clipped && vertex.isClipped()) {
      // CASE WHERE STARTING VERTEX IS CLIPPED
      auto current_edge = vertex.getHalfEdge();
      const auto starting_edge = current_edge;
      do {
        // If it has needsToSeek set, it means it is a new vertex
        // or already visited. Either way, do not need to check
        // for intersection
        const auto& vertex_start = current_edge->getPreviousVertex();
        if (vertex_start->needsToSeek()) {
          current_edge =
              current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
          continue;
        }
        const auto& vertex_end = current_edge->getVertex();
        if (elliptic) {
          if (a_aligned_paraboloid.a() > ZERO) {
            if (vertex_start->isNotClipped() && vertex_end->isNotClipped()) {
              current_edge =
                  current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
              continue;
            }
          } else if (vertex_start->isClipped() && vertex_end->isClipped()) {
            current_edge =
                current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
            continue;
          }
        }

        // If previous vertex is not clipped, single-intercept
        // If previous vertex is clipped, need to check for
        // double-intercept Checking for double-intercept and
        // calculating single-intercept is the same routine, so
        // just always do it.
        const auto& edge_start = vertex_start->getLocation();
        const auto& edge_end = vertex_end->getLocation();
        checkAndFindIntercepts<pt_type, ScalarType>(
            a_aligned_paraboloid, edge_start, edge_end, &edge_intercepts,
            nudge_epsilon, elliptic);

        // assert(current_edge->getPreviousVertex()->isClipped()
        // ||
        //        edge_intercepts.size() == 1);
        // Size of returned intercepts indicates single or double
        // intercept (or none)
        if (edge_intercepts.size() == 1) {
          // Check for intersection near end point
          if (squaredMagnitude(edge_intercepts[0] - edge_start) <
                  nudge_epsilon_sq ||
              squaredMagnitude(edge_intercepts[0] - edge_end) <
                  nudge_epsilon_sq) {
            requires_nudge = true;
            break;
          }

          // Add vertex to half edge structure, resetting
          // connectivity and creating a new half edge (and new
          // opposite half edge)
          placeSingleIntercept(edge_intercepts[0], current_edge, a_polytope,
                               a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->addIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->setStartingHalfEdge(opposite_half_edge);
          opposite_face->addIntersection();
        } else if (edge_intercepts.size() == 2) {
          // Check for intersection near end point
          if (squaredMagnitude(edge_intercepts[0] - edge_start) <
                  nudge_epsilon_sq ||
              squaredMagnitude(edge_intercepts[1] - edge_end) <
                  nudge_epsilon_sq ||
              squaredMagnitude(edge_intercepts[0] - edge_intercepts[1]) <
                  nudge_epsilon_sq) {
            requires_nudge = true;
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
        }
        current_edge =
            current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
      } while (current_edge != starting_edge);
    } else if (check_from_unclipped && vertex.isNotClipped()) {
      // CASE WHERE STARTING VERTEX IS UNCLIPPED
      auto current_edge = vertex.getHalfEdge();
      const auto starting_edge = current_edge;
      do {
        // If it has needsToSeek set, it means it is a new vertex
        // or already visited. Either way, do not need to check
        // for intersection
        const auto& vertex_start = current_edge->getPreviousVertex();
        if (vertex_start->needsToSeek()) {
          current_edge =
              current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
          continue;
        }
        const auto& vertex_end = current_edge->getVertex();
        if (elliptic) {
          if (a_aligned_paraboloid.a() > ZERO) {
            if (vertex_start->isNotClipped() && vertex_end->isNotClipped()) {
              current_edge =
                  current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
              continue;
            }
          } else if (vertex_start->isClipped() && vertex_end->isClipped()) {
            current_edge =
                current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
            continue;
          }
        }

        // If previous vertex is clipped, single-intercept
        // If previous vertex is not clipped, need to check for
        // double-intercept Checking for double-intercept and
        // calculating single-intercept is the same routine, so
        // just always do it.
        const auto& edge_start = vertex_start->getLocation();
        const auto& edge_end = vertex_end->getLocation();
        checkAndFindIntercepts<pt_type, ScalarType>(
            a_aligned_paraboloid, edge_start, edge_end, &edge_intercepts,
            nudge_epsilon, elliptic);

        // assert(current_edge->getPreviousVertex()->isNotClipped()
        // ||
        //        edge_intercepts.size() == 1);
        // Size of returned intercepts indicates single or double
        // intercept (or none)
        if (edge_intercepts.size() == 1) {
          // Check for intersection near end point
          if (squaredMagnitude(edge_intercepts[0] - edge_start) <
                  nudge_epsilon_sq ||
              squaredMagnitude(edge_intercepts[0] - edge_end) <
                  nudge_epsilon_sq) {
            requires_nudge = true;
            break;
          }

          // Add vertex to half edge structure, resetting
          // connectivity and creating a new half edge (and new
          // opposite half edge)
          placeSingleIntercept(edge_intercepts[0], current_edge, a_polytope,
                               a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->setStartingHalfEdge(
              current_edge->getPreviousHalfEdge());
          current_face->addIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->addIntersection();
        } else if (edge_intercepts.size() == 2) {
          // Check for intersection near end point
          if (squaredMagnitude(edge_intercepts[0] - edge_start) <
                  nudge_epsilon_sq ||
              squaredMagnitude(edge_intercepts[1] - edge_end) <
                  nudge_epsilon_sq ||
              squaredMagnitude(edge_intercepts[0] - edge_intercepts[1]) <
                  nudge_epsilon_sq) {
            requires_nudge = true;
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
  // After leaving above loop, all faces that were intersected
  // have a starting half edge that ends on a new (intersection)
  // vertex and the previous vertex is unclipped (or new). (i.e.,
  // the edge exists in the unclipped portion. All new vertices
  // also have their half edge being one that has its edge in the
  // unclipped portion.

  const auto vertices_after_intersection = a_polytope->getNumberOfVertices();
  const auto new_intersection_vertices =
      vertices_after_intersection - starting_number_of_vertices;

  for (UnsignedIndex_t v = starting_number_of_vertices;
       v < vertices_after_intersection; ++v) {
    // Original vertices will be set as needsToSeek()
    // Reset newly created vertices will have doesNotNeedToSeek()
    a_polytope->getVertex(v)->setAsUnnecessaryToSeek();
  }

  /* Triangulate faces with intersections */
  if (elliptic) {
    triangulatePolytope(a_polytope, a_complete_polytope, a_aligned_paraboloid,
                        nudge_epsilon, &requires_nudge);
  }

  // If intersection is too close to corners, shift polyhedron and
  // try again
  if (requires_nudge) {
    // Clean half-edge structure by removing intersections
    assert(a_polytope->checkValidHalfEdgeStructure());
    if (new_intersection_vertices > 0) {
      resetPolyhedron(a_polytope, a_complete_polytope);
    }
    assert(a_polytope->checkValidHalfEdgeStructure());

    AlignedParaboloidBase<Quad_t> nudged_aligned_paraboloid =
        nudgeParaboloid(a_aligned_paraboloid, a_nudge_iter + 1);

    if constexpr (!has_embedded_gradient<ReturnType>::value) {
      if constexpr (std::is_same_v<ScalarType, double>) {
        converted_complete_polytope_type converted_polytope_paraboloid;
        convertPolytopeFromDoubleToQuadPrecision(
            a_polytope, a_complete_polytope, &converted_polytope_paraboloid);
        auto converted_segmented_paraboloid =
            converted_polytope_paraboloid.generateSegmentedPolyhedron();
        assert(converted_segmented_paraboloid.checkValidHalfEdgeStructure());

        /* Nudge polytope */
        nudgePolyhedron(&converted_segmented_paraboloid,
                        &converted_polytope_paraboloid, a_nudge_iter + 1,
                        a_surface);
        return formParaboloidIntersectionBases<ReturnType>(
            &converted_segmented_paraboloid, &converted_polytope_paraboloid,
            nudged_aligned_paraboloid, a_nudge_iter + 1, a_surface);
      } else {
        nudgePolyhedron(a_polytope, a_complete_polytope, a_nudge_iter + 1,
                        a_surface);
        return formParaboloidIntersectionBases<ReturnType>(
            a_polytope, a_complete_polytope, nudged_aligned_paraboloid,
            a_nudge_iter + 1, a_surface);
      }
    } else {
      ///////////////////////////////////////TODO }
    }
  }

  // Recalculate in face plane information
  for (auto& face : (*a_polytope)) {
    auto current_edge = getNextNonIntersection(face->getStartingHalfEdge());
    const auto& p0 = current_edge->getVertex()->getLocation().getPt();
    current_edge = getNextNonIntersection(current_edge->getNextHalfEdge());
    const auto& p1 = current_edge->getVertex()->getLocation().getPt();
    current_edge = getNextNonIntersection(current_edge->getNextHalfEdge());
    const auto& p2 = current_edge->getVertex()->getLocation().getPt();
    auto normal = Normal(
        (p1[1] - p0[1]) * (p2[2] - p1[2]) - (p1[2] - p0[2]) * (p2[1] - p1[1]),
        (p1[2] - p0[2]) * (p2[0] - p1[0]) - (p1[0] - p0[0]) * (p2[2] - p1[2]),
        (p1[0] - p0[0]) * (p2[1] - p1[1]) - (p1[1] - p0[1]) * (p2[0] - p1[0]));
    ScalarType squaredMag = squaredMagnitude(normal);
    if (squaredMag < DISTANCE_EPSILON * DISTANCE_EPSILON) {
      normal = Normal(ZERO, ZERO, ZERO);
    } else {
      normal /= sqrt(squaredMag);
    }
    face->setPlane(Plane(normal, normal * p0));
  }

  // Check for face-only intersections. Can only happen with
  // elliptic
  if (elliptic) {
    if (a_aligned_paraboloid.a() > ZERO) {
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
          if (vertex.getLocation().getPt()[2] < ZERO) {
            face_valid = true;
          }
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != starting_half_edge);
        if (!face_valid) {
          continue;
        }

        // Made it this far, check if there is a face-only
        // intersection
        const auto& face_plane = face.getPlane();
        if (fabs(face_plane.normal()[2]) > MACHINE_EPSILON) {
          // Get ellipse on this face
          if (ellipseContainedInFace<ScalarType>(
                  a_aligned_paraboloid, face_plane, starting_half_edge)) {
            if constexpr (!has_embedded_gradient<ReturnType>::value) {
              full_moments +=
                  computeFaceOnlyContribution<ReturnType, ScalarType>(
                      a_aligned_paraboloid, face_plane,
                      starting_half_edge->getVertex()->getLocation());
            } else {
              full_moments +=
                  computeFaceOnlyContributionWithGradient<ReturnType,
                                                          ScalarType>(
                      a_aligned_paraboloid, face_plane,
                      starting_half_edge->getVertex()->getLocation());
            }
            // Return surface parametrization
            if constexpr (!std::is_same<SurfaceOutputType,
                                        NoSurfaceOutput>::value) {
              addEllipseToSurfaceOutput<ScalarType>(a_aligned_paraboloid,
                                                    face_plane, a_surface);
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
          if (vertex.getLocation().getPt()[2] > ZERO) {
            face_valid = true;
          }
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != starting_half_edge);
        if (!face_valid) {
          continue;
        }

        // Made it this far, check if there is a face-only
        // intersection
        const auto& face_plane = face.getPlane();
        if (fabs(face_plane.normal()[2]) > MACHINE_EPSILON) {
          // Get ellipse on this face
          if (ellipseContainedInFace<ScalarType>(
                  a_aligned_paraboloid, face_plane, starting_half_edge)) {
            if constexpr (!has_embedded_gradient<ReturnType>::value) {
              full_moments +=
                  computeFaceOnlyContribution<ReturnType, ScalarType>(
                      a_aligned_paraboloid, face_plane,
                      starting_half_edge->getVertex()->getLocation());
            } else {
              full_moments +=
                  computeFaceOnlyContributionWithGradient<ReturnType,
                                                          ScalarType>(
                      a_aligned_paraboloid, face_plane,
                      starting_half_edge->getVertex()->getLocation());
            }
            if constexpr (!std::is_same<SurfaceOutputType,
                                        NoSurfaceOutput>::value) {
              addEllipseToSurfaceOutput<ScalarType>(a_aligned_paraboloid,
                                                    face_plane, a_surface);
            }
          }
        }
      }
    }
  }

  if (new_intersection_vertices == 0) {
    // If no intersections at all, gradient need to be updated
    if constexpr (has_embedded_gradient<ReturnType>::value) {
      // if constexpr (std::is_same<typename
      // ReturnType::moment_type,
      //                            VolumeMoments>::value) {
      // If volume == 0, then there are no face-only intersections
      if (full_moments.volume() == ZERO) {
        // Find vertex closest to the paraboloid
        ScalarType min_dist = static_cast<ScalarType>(DBL_MAX);
        pt_type closest_pt;
        for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
          auto& vertex = *(a_polytope->getVertex(v));
          const ScalarType dist = fabs(signedDistance(
              vertex.getLocation().getPt(), a_aligned_paraboloid));
          if (dist < min_dist) {
            min_dist = dist;
            closest_pt = vertex.getLocation();
          }
        }
        // Add (normalized) signed distance to the volume
        using gradient_type = typename ReturnType::gradient_type;
        using scalar_type = ScalarWithGradient<gradient_type>;
        std::array<scalar_type, 3> pt;
        for (UnsignedIndex_t d = 0; d < 3; ++d) {
          pt[d] = scalar_type(closest_pt.getPt()[d], closest_pt.getData()[d]);
        }
        auto A_grad = gradient_type(ZERO), B_grad = gradient_type(ZERO);
        A_grad.setGradA(ONE);
        B_grad.setGradB(ONE);
        const auto A = scalar_type(a_aligned_paraboloid.a(), A_grad),
                   B = scalar_type(a_aligned_paraboloid.b(), B_grad);
        const auto volume_add_on =
            -(A * pt[0] * pt[0] + B * pt[1] * pt[1] + pt[2]);
        const ScalarType scale =
            pow(static_cast<ScalarType>(a_polytope->calculateVolume()),
                TWO / THREE);
        full_moments.volume() = volume_add_on.value() * scale;
        full_moments.volume_gradient() = volume_add_on.gradient() * scale;
      }
    }

    if (number_of_vertices_above == starting_number_of_vertices) {
      // All points above
      return full_moments;
    } else {
      // All points below
      return ReturnType::calculateMoments(a_polytope) + full_moments;
    }
  }

  for (UnsignedIndex_t f = 0; f < a_polytope->getNumberOfFaces(); ++f) {
    auto& face = *(*a_polytope)[f];
    const auto& face_normal = face.getPlane().normal();
    if (squaredMagnitude(face_normal) < MACHINE_EPSILON) {
      continue;
    }
    auto starting_half_edge = face.getStartingHalfEdge();
    const auto intersection_size = face.getNumberOfIntersections();
    if (intersection_size % 2 == 1) {
      resetPolyhedron(a_polytope, a_complete_polytope);

      AlignedParaboloidBase<Quad_t> nudged_aligned_paraboloid =
          nudgeParaboloid(a_aligned_paraboloid, a_nudge_iter + 1);

      if constexpr (!has_embedded_gradient<ReturnType>::value) {
        if constexpr (std::is_same_v<ScalarType, double>) {
          converted_complete_polytope_type converted_polytope_paraboloid;
          convertPolytopeFromDoubleToQuadPrecision(
              a_polytope, a_complete_polytope, &converted_polytope_paraboloid);
          auto converted_segmented_paraboloid =
              converted_polytope_paraboloid.generateSegmentedPolyhedron();
          assert(converted_segmented_paraboloid.checkValidHalfEdgeStructure());

          /* Nudge polytope */
          nudgePolyhedron(&converted_segmented_paraboloid,
                          &converted_polytope_paraboloid, a_nudge_iter + 1,
                          a_surface);
          return formParaboloidIntersectionBases<ReturnType>(
              &converted_segmented_paraboloid, &converted_polytope_paraboloid,
              nudged_aligned_paraboloid, a_nudge_iter + 1, a_surface);
        } else {
          nudgePolyhedron(a_polytope, a_complete_polytope, a_nudge_iter + 1,
                          a_surface);
          return formParaboloidIntersectionBases<ReturnType>(
              a_polytope, a_complete_polytope, nudged_aligned_paraboloid,
              a_nudge_iter + 1, a_surface);
        }
      } else {
        ///////////////////////////////////////TODO }
      }
    }

    if (intersection_size == 0) {
      if (starting_half_edge->getVertex()->isNotClipped()) {
        // The face is entirely below
        const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
        auto current_half_edge =
            starting_half_edge->getNextHalfEdge()->getNextHalfEdge();
        auto prev_pt = current_half_edge->getPreviousVertex()->getLocation();
        do {
          const auto& curr_pt = current_half_edge->getVertex()->getLocation();
          if constexpr (!has_embedded_gradient<ReturnType>::value) {
            full_moments += computeType1Contribution<ReturnType, ScalarType>(
                ref_pt, prev_pt, curr_pt);
          } else {
            full_moments +=
                computeType1ContributionWithGradient<ReturnType, ScalarType>(
                    ref_pt, prev_pt, curr_pt);
          }
          prev_pt = curr_pt;
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != starting_half_edge);
      }
    } else if (intersection_size == 2) {
      const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
      half_edge_type* exit_half_edge;
      full_moments +=
          computeUnclippedSegmentType1Contribution<ReturnType, ScalarType>(
              a_aligned_paraboloid, ref_pt, starting_half_edge, exit_half_edge,
              true);
      full_moments += computeNewEdgeSegmentContribution<ReturnType, ScalarType>(
          a_aligned_paraboloid, ref_pt, starting_half_edge, exit_half_edge,
          true, false, &requires_nudge, a_surface);
    } else {
      // More than 2 intersections. More complicated cases.
      const bool elliptic_face =
          a_aligned_paraboloid.a() * a_aligned_paraboloid.b() > ZERO &&
          fabs(face_normal[2]) > MACHINE_EPSILON;
      const bool hyperbolic_face =
          a_aligned_paraboloid.a() * a_aligned_paraboloid.b() < ZERO &&
          fabs(face_normal[2]) > MACHINE_EPSILON;

      // -----> NO SORTING
      // if constexpr (std::is_same<SurfaceOutputType, NoSurfaceOutput>::value)
      // {
      if (elliptic || std::is_same_v<SurfaceOutputType, NoSurfaceOutput>) {
        if (elliptic_face) {
          const bool reverse = a_aligned_paraboloid.a() < ZERO;
          const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
          half_edge_type* exit_half_edge;
          auto current_edge = starting_half_edge;
          bool skip_first = true;
          std::size_t found_intersections = 0;
          do {
            full_moments +=
                computeUnclippedSegmentType1Contribution<ReturnType,
                                                         ScalarType>(
                    a_aligned_paraboloid, ref_pt, current_edge, exit_half_edge,
                    skip_first);

            if (reverse) {
              full_moments +=
                  computeNewEdgeSegmentContribution<ReturnType, ScalarType>(
                      a_aligned_paraboloid, ref_pt, current_edge,
                      exit_half_edge, false, false, &requires_nudge, a_surface);
              current_edge = exit_half_edge->getNextHalfEdge();
              while (current_edge->getVertex()->needsToSeek()) {
                current_edge = current_edge->getNextHalfEdge();
              }
            } else {
              current_edge = exit_half_edge->getNextHalfEdge();
              while (current_edge->getVertex()->needsToSeek()) {
                current_edge = current_edge->getNextHalfEdge();
              }
              full_moments +=
                  computeNewEdgeSegmentContribution<ReturnType, ScalarType>(
                      a_aligned_paraboloid, ref_pt, current_edge,
                      exit_half_edge, false, false, &requires_nudge, a_surface);
            }
            skip_first = false;
            found_intersections += 2;
          } while (found_intersections != intersection_size);
        } else {
          SmallVector<half_edge_type*, 4> intersections;
          // Find intersections and determine status
          auto current_edge = starting_half_edge;
          bool reverse_order = false;
          half_edge_type* exit_half_edge;
          std::size_t found_intersections = 0;
          const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
          bool skip_first = true;
          do {
            full_moments +=
                computeUnclippedSegmentType1Contribution<ReturnType,
                                                         ScalarType>(
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

          // if (elliptic_face) {
          //   reverse_order = a_aligned_paraboloid.a() < ZERO;
          // } else
          if (hyperbolic_face) {
            // std::cout << "HYPERBOLIC FACE with more than 2
            // intersections"
            // << std::endl;
            const std::array<ScalarType, 2> conic_center{
                {face_normal[0] /
                     (TWO * a_aligned_paraboloid.a() * face_normal[2]),
                 face_normal[1] /
                     (TWO * a_aligned_paraboloid.b() * face_normal[2])}};
            const auto& pt_in_plane =
                intersections[0]->getVertex()->getLocation().getPt();
            const ScalarType delta_face = (face_normal[0] * pt_in_plane[0] +
                                           face_normal[1] * pt_in_plane[1] +
                                           face_normal[2] * pt_in_plane[2]) /
                                          face_normal[2];
            const ScalarType gamma_face =
                a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
                a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] -
                delta_face;
            const ScalarType z_diff =
                face_normal[0] * conic_center[0] / face_normal[2] +
                face_normal[1] * conic_center[1] / face_normal[2] - gamma_face -
                TWO * delta_face;
            const std::size_t split_ind =
                a_aligned_paraboloid.a() * gamma_face > ZERO ? 0 : 1;

            bool vertices_on_same_branch = true;
            for (std::size_t i = 1; i < intersection_size; ++i) {
              const auto& curr_pt =
                  intersections[i]->getVertex()->getLocation().getPt();
              if ((pt_in_plane[split_ind] - conic_center[split_ind]) *
                      (curr_pt[split_ind] - conic_center[split_ind]) <
                  ZERO) {
                vertices_on_same_branch = false;
                break;
              }
            }
            bool total_invert = z_diff < ZERO ? !vertices_on_same_branch
                                              : vertices_on_same_branch;
            if (total_invert) {
              reverse_order = true;
            }
          } else {
            // std::cout << "PARABOLIC FACE with more than 2
            // intersections"
            //           << std::endl;

            Pt intersection_avg = Pt(ZERO, ZERO, ZERO);
            for (std::size_t i = 0; i < intersection_size; ++i) {
              const auto& curr_pt =
                  intersections[i]->getVertex()->getLocation().getPt();
              intersection_avg += curr_pt;
            }
            intersection_avg /= static_cast<ScalarType>(intersection_size);
            const ScalarType z_diff =
                intersection_avg[2] +
                a_aligned_paraboloid.a() * intersection_avg[0] *
                    intersection_avg[0] +
                a_aligned_paraboloid.b() * intersection_avg[1] *
                    intersection_avg[1];
            if (z_diff > ZERO) {
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
            full_moments +=
                computeNewEdgeSegmentContribution<ReturnType, ScalarType>(
                    a_aligned_paraboloid, ref_pt, entry_half_edge,
                    exit_half_edge, first_entry, false, &requires_nudge,
                    a_surface);
            first_entry = false;
          }
          // Clear list of intersections
          intersections.clear();
        }
      }
      // -- --->If surface ouput : NEEDS SORTING
      else {
        bool ignore_type3_contributions = false;
        using stype = std::pair<half_edge_type*, ScalarType>;
        SmallVector<stype, 10> intersections;
        // Find intersections and determine status
        auto current_edge = starting_half_edge;
        bool reverse_order = false;
        half_edge_type* exit_half_edge;
        std::size_t found_intersections = 0;
        const auto& ref_pt = starting_half_edge->getVertex()->getLocation();
        bool skip_first = true;
        do {
          full_moments +=
              computeUnclippedSegmentType1Contribution<ReturnType, ScalarType>(
                  a_aligned_paraboloid, ref_pt, current_edge, exit_half_edge,
                  skip_first);
          current_edge->getVertex()->markAsEntry();
          exit_half_edge->getVertex()->markAsExit();
          intersections.push_back(std::pair<half_edge_type*, ScalarType>(
              {current_edge, static_cast<ScalarType>(DBL_MAX)}));
          intersections.push_back(std::pair<half_edge_type*, ScalarType>(
              {exit_half_edge, static_cast<ScalarType>(DBL_MAX)}));
          current_edge = exit_half_edge->getNextHalfEdge();
          while (current_edge->getVertex()->needsToSeek()) {
            current_edge = current_edge->getNextHalfEdge();
          }
          found_intersections += 2;
          skip_first = false;
        } while (found_intersections != intersection_size);

        if (elliptic_face) {
          const std::array<ScalarType, 2> conic_center{
              {face_normal[0] /
                   (TWO * a_aligned_paraboloid.a() * face_normal[2]),
               face_normal[1] /
                   (TWO * a_aligned_paraboloid.b() * face_normal[2])}};
          const ScalarType normal_invert = copysign(ONE, face_normal[2]);
          const ScalarType invert =
              a_aligned_paraboloid.a() < ZERO ? -normal_invert : normal_invert;

          for (auto& element : intersections) {
            const auto& pt = element.first->getVertex()->getLocation().getPt();
            element.second = invert * atan2(pt[1] - conic_center[1],
                                            pt[0] - conic_center[0]);
          }
          std::sort(intersections.begin(), intersections.end(),
                    [](const stype& a, const stype& b) {
                      return a.second < b.second;
                    });

          // Backdoor for case where two angles are equal (because
          // of floating-point errors)
          bool restart_sort = false;
          for (std::size_t i = 0; i < intersection_size - 1; ++i) {
            if (fabs(intersections[i].second - intersections[i + 1].second) <
                ONE_HUNDRED * MACHINE_EPSILON) {
              restart_sort = true;
              break;
            }
          }
          restart_sort = true;
          if (restart_sort) {
            auto intersection_copy = intersections;
            // We split in X direction
            UnsignedIndex_t split_ind = 0;
            UnsignedIndex_t store_ind = 1;
            UnsignedIndex_t pos_end = 0;
            UnsignedIndex_t neg_end = 0;
            for (auto& element : intersections) {
              const auto& pt =
                  element.first->getVertex()->getLocation().getPt();
              if (pt[split_ind] > conic_center[split_ind]) {
                const auto loc = pos_end++;
                intersection_copy[loc].first = element.first;
                intersection_copy[loc].second = invert * pt[store_ind];
              } else {
                const auto loc = intersection_size - 1 - (neg_end++);
                intersection_copy[loc].first = element.first;
                intersection_copy[loc].second = -invert * pt[store_ind];
              }
            }
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
            bool re_restart_sort = false;
            if (pos_end > 0) {
              for (UnsignedIndex_t i = 0; i < pos_end - 1; ++i) {
                if (fabs(intersections[i].second -
                         intersections[i + 1].second) <
                    ONE_HUNDRED * MACHINE_EPSILON) {
                  re_restart_sort = true;
                  break;
                }
              }
            }
            for (UnsignedIndex_t i = pos_end; i < intersection_size - 1; ++i) {
              if (fabs(intersections[i].second - intersections[i + 1].second) <
                  ONE_HUNDRED * MACHINE_EPSILON) {
                re_restart_sort = true;
                break;
              }
            }
            if (re_restart_sort) {
              // We split in Y direction
              split_ind = 1;
              store_ind = 0;
              pos_end = 0;
              neg_end = 0;
              for (auto& element : intersections) {
                const auto& pt =
                    element.first->getVertex()->getLocation().getPt();
                if (pt[split_ind] > conic_center[split_ind]) {
                  const auto loc = pos_end++;
                  intersection_copy[loc].first = element.first;
                  intersection_copy[loc].second = -invert * pt[store_ind];
                } else {
                  const auto loc = intersection_size - 1 - (neg_end++);
                  intersection_copy[loc].first = element.first;
                  intersection_copy[loc].second = invert * pt[store_ind];
                }
              }
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
              bool re_re_restart_sort = false;
              if (pos_end > 0) {
                for (UnsignedIndex_t i = 0; i < pos_end - 1; ++i) {
                  if (fabs(intersections[i].second -
                           intersections[i + 1].second) <
                      ONE_HUNDRED * MACHINE_EPSILON) {
                    re_re_restart_sort = true;
                    break;
                  }
                }
              }
              for (UnsignedIndex_t i = pos_end; i < intersection_size - 1;
                   ++i) {
                if (fabs(intersections[i].second -
                         intersections[i + 1].second) <
                    ONE_HUNDRED * MACHINE_EPSILON) {
                  re_re_restart_sort = true;
                  break;
                }
              }
              if (re_re_restart_sort) {
                resetPolyhedron(a_polytope, a_complete_polytope);

                if constexpr (!has_embedded_gradient<ReturnType>::value) {
                  if constexpr (std::is_same_v<ScalarType, double>) {
                    converted_complete_polytope_type
                        converted_polytope_paraboloid;
                    convertPolytopeFromDoubleToQuadPrecision(
                        a_polytope, a_complete_polytope,
                        &converted_polytope_paraboloid);
                    auto converted_segmented_paraboloid =
                        converted_polytope_paraboloid
                            .generateSegmentedPolyhedron();
                    assert(converted_segmented_paraboloid
                               .checkValidHalfEdgeStructure());

                    /* Nudge polytope */
                    nudgePolyhedron(&converted_segmented_paraboloid,
                                    &converted_polytope_paraboloid,
                                    a_nudge_iter + 1, a_surface);
                    return formParaboloidIntersectionBases<ReturnType>(
                        &converted_segmented_paraboloid,
                        &converted_polytope_paraboloid,
                        AlignedParaboloidBase<Quad_t>(a_aligned_paraboloid),
                        a_nudge_iter + 1, a_surface);
                  } else {
                    nudgePolyhedron(a_polytope, a_complete_polytope,
                                    a_nudge_iter + 1, a_surface);
                    return formParaboloidIntersectionBases<ReturnType>(
                        a_polytope, a_complete_polytope, a_aligned_paraboloid,
                        a_nudge_iter + 1, a_surface);
                  }
                } else {
                  ///////////////////////////////////////TODO }
                }
              }
            }
          }
        } else if (hyperbolic_face) {
          std::size_t pos_end = 0;
          std::size_t neg_end = 0;
          auto intersection_copy = intersections;
          const std::array<ScalarType, 2> conic_center{
              {face_normal[0] /
                   (TWO * a_aligned_paraboloid.a() * face_normal[2]),
               face_normal[1] /
                   (TWO * a_aligned_paraboloid.b() * face_normal[2])}};
          const auto& pt_in_plane =
              intersections[0].first->getVertex()->getLocation().getPt();
          const ScalarType delta_face = (face_normal[0] * pt_in_plane[0] +
                                         face_normal[1] * pt_in_plane[1] +
                                         face_normal[2] * pt_in_plane[2]) /
                                        face_normal[2];
          const ScalarType gamma_face =
              a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
              a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] -
              delta_face;
          const std::size_t split_ind =
              a_aligned_paraboloid.a() * gamma_face > ZERO ? 0 : 1;
          const std::size_t store_ind = split_ind == 0 ? 1 : 0;
          const ScalarType z_center_plane =
              -face_normal[0] * conic_center[0] / face_normal[2] -
              face_normal[1] * conic_center[1] / face_normal[2] + delta_face;
          const ScalarType z_center_paraboloid =
              -a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] -
              a_aligned_paraboloid.b() * conic_center[1] * conic_center[1];
          ScalarType total_invert = copysign(ONE, face_normal[2]);
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
          if (fabs(gamma_face) < ONE_HUNDRED * MACHINE_EPSILON ||
              fabs(z_center_plane - z_center_paraboloid) <
                  ONE_HUNDRED * MACHINE_EPSILON) {
            requires_nudge = true;
          } else {
            for (std::size_t i = 0; i < intersection_size - 1; ++i) {
              if (fabs(intersection_copy[i].second -
                       intersection_copy[i + 1].second) <
                  ONE_HUNDRED * MACHINE_EPSILON) {
                requires_nudge = true;
                break;
              }
            }
          }
          intersections = intersection_copy;
        } else {
          Pt intersection_avg = Pt(ZERO, ZERO, ZERO);
          for (std::size_t i = 0; i < intersection_size; ++i) {
            const auto& curr_pt =
                intersections[i].first->getVertex()->getLocation().getPt();
            intersection_avg += curr_pt;
          }
          intersection_avg /= static_cast<ScalarType>(intersection_size);
          const ScalarType z_diff =
              intersection_avg[2] +
              a_aligned_paraboloid.a() * intersection_avg[0] *
                  intersection_avg[0] +
              a_aligned_paraboloid.b() * intersection_avg[1] *
                  intersection_avg[1];

          if (fabs(z_diff) < ONE_HUNDRED * MACHINE_EPSILON) {
            requires_nudge = true;
          }

          // Find dominant face normal direction
          std::size_t dir = 0;
          ScalarType max_normal = fabs(face_normal[dir]);
          for (std::size_t d = 1; d < 3; ++d) {
            if (fabs(face_normal[d]) > max_normal) {
              dir = d;
            }
            max_normal = fabs(face_normal[dir]);
          }
          const ScalarType normal_invert = copysign(ONE, face_normal[dir]);
          const ScalarType invert =
              z_diff < ZERO ? normal_invert : -normal_invert;
          // Compute convex hull on projected plane
          const std::size_t x_id = (dir + 1) % 3;
          const std::size_t y_id = (dir + 2) % 3;

          // Find the leftmost point
          std::size_t left_id = 0;
          auto left_pt =
              intersections[left_id].first->getVertex()->getLocation().getPt();
          for (UnsignedIndex_t i = 1; i < intersection_size; ++i) {
            const auto pt =
                intersections[i].first->getVertex()->getLocation().getPt();
            if (pt[x_id] < left_pt[x_id]) {
              left_id = i;
              left_pt = pt;
            }
          }

          // Start from leftmost point, keep moving
          // counterclockwise until reach the start point again.
          std::size_t p = left_id;
          auto p_pt = left_pt;
          std::size_t hull_size = 0;
          auto intersection_copy = intersections;
          bool is_flat = false;
          bool has_flat = false;
          do {
            p_pt = intersections[p].first->getVertex()->getLocation().getPt();
            intersection_copy[hull_size++] = intersections[p];
            std::size_t q = (p + 1) % intersection_size;
            std::size_t flatness_counter = 0;
            for (std::size_t i = 0; i < intersection_size; ++i) {
              if (q != i && p != i) {
                const auto i_pt =
                    intersections[i].first->getVertex()->getLocation().getPt();
                const auto q_pt =
                    intersections[q].first->getVertex()->getLocation().getPt();
                ScalarType dot_product =
                    (i_pt[y_id] - p_pt[y_id]) * (q_pt[x_id] - i_pt[x_id]) -
                    (i_pt[x_id] - p_pt[x_id]) * (q_pt[y_id] - i_pt[y_id]);
                // If the points are aligned, keep closest and
                // increase flatness counter
                if (fabs(dot_product) < ANGLE_EPSILON) {
                  has_flat = true;
                  flatness_counter++;
                } else if (dot_product < ZERO) {
                  q = i;
                }
              }
            }
            // Check if all intersections are aligned
            if (flatness_counter == intersection_size - 2) {
              is_flat = true;
            }
            p = q;
          } while (p != left_id && !is_flat && hull_size <= intersection_size);
          if (!is_flat) {
            if (!is_flat && has_flat) {
              for (auto& element : intersections) {
                const auto& pt =
                    element.first->getVertex()->getLocation().getPt();
                // element.second =
                //     invert *
                //     std::atan2(pt[store_ind] -
                //     intersection_avg[store_ind],
                //                pt[x_id] -
                //                intersection_avg[x_id]);
                element.second =
                    invert *
                    copysign(
                        ONE - (pt[x_id] - intersection_avg[x_id]) /
                                  (fabs(pt[x_id] - intersection_avg[x_id]) +
                                   fabs(pt[y_id] - intersection_avg[y_id])),
                        (pt[y_id] - intersection_avg[y_id]));
              }
              std::sort(intersections.begin(), intersections.end(),
                        [](const stype& a, const stype& b) {
                          return a.second < b.second;
                        });
            } else {
              if (hull_size != intersection_size) {
                std::cout << "Parabolic intersections don't form a "
                             "convex polygon!"
                          << std::endl;
                exit(-1);
              }
              if (invert > ZERO) {
                intersections = intersection_copy;
              } else {
                for (std::size_t i = 0; i < intersection_size; ++i) {
                  intersections[i] =
                      intersection_copy[intersection_size - 1 - i];
                }
              }
            }
          }
          UnsignedIndex_t intersection_are_aligned = 0;
          for (std::size_t i = 0; i < intersection_size - 1; ++i) {
            if (fabs(intersections[i].second - intersections[i + 1].second) <
                ONE_HUNDRED * MACHINE_EPSILON) {
              requires_nudge = true;
              intersection_are_aligned++;  // break;
            }
          }
          if (intersection_are_aligned == intersection_size) {
            ignore_type3_contributions = true;
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
          full_moments +=
              computeNewEdgeSegmentContribution<ReturnType, ScalarType>(
                  a_aligned_paraboloid, ref_pt, entry_half_edge, exit_half_edge,
                  false, ignore_type3_contributions, &requires_nudge,
                  a_surface);
        }
        // Clear list of intersections
        intersections.clear();
      }
    }
    if (requires_nudge) {
      resetPolyhedron(a_polytope, a_complete_polytope);

      AlignedParaboloidBase<Quad_t> nudged_aligned_paraboloid =
          nudgeParaboloid(a_aligned_paraboloid, a_nudge_iter + 1);

      if constexpr (!has_embedded_gradient<ReturnType>::value) {
        if constexpr (std::is_same_v<ScalarType, double>) {
          converted_complete_polytope_type converted_polytope_paraboloid;
          convertPolytopeFromDoubleToQuadPrecision(
              a_polytope, a_complete_polytope, &converted_polytope_paraboloid);
          auto converted_segmented_paraboloid =
              converted_polytope_paraboloid.generateSegmentedPolyhedron();
          assert(converted_segmented_paraboloid.checkValidHalfEdgeStructure());

          /* Nudge polytope */
          nudgePolyhedron(&converted_segmented_paraboloid,
                          &converted_polytope_paraboloid, a_nudge_iter + 1,
                          a_surface);
          return formParaboloidIntersectionBases<ReturnType>(
              &converted_segmented_paraboloid, &converted_polytope_paraboloid,
              nudged_aligned_paraboloid, a_nudge_iter + 1, a_surface);
        } else {
          nudgePolyhedron(a_polytope, a_complete_polytope, a_nudge_iter + 1,
                          a_surface);
          return formParaboloidIntersectionBases<ReturnType>(
              a_polytope, a_complete_polytope, nudged_aligned_paraboloid,
              a_nudge_iter + 1, a_surface);
        }
      } else {
        ///////////////////////////////////////TODO }
      }
    }
  }

  // //////////////////////// FOR FIGURE
  // if constexpr (!std::is_same<SurfaceOutputType,
  // NoSurfaceOutput>::value)
  // {
  //   std::cout << std::endl
  //             << std::endl
  //             << "/////// START ASYMPTOTE CODE" << std::endl;
  //   std::cout << "real scale_vertex = 0.02;" << std::endl;
  //   std::cout << "real scale_edge = 0.005;" << std::endl;
  //   std::cout << "real scale_arc = 0.006;" << std::endl;
  //   std::cout << "pen color_corner_unclipped = black;" <<
  //   std::endl; std::cout << "pen color_proj_corner_unclipped =
  //   black;" << std::endl; std::cout << "pen
  //   color_corner_clipped = gray;" << std::endl; std::cout <<
  //   "pen color_intersection = red;" << std::endl; std::cout <<
  //   "pen color_proj_intersection = red;" << std::endl;
  //   std::cout << "pen color_edge = gray;" << std::endl;
  //   std::cout << "pen color_unclipped_edge = black;" <<
  //   std::endl; std::cout << "pen color_arc = red;" <<
  //   std::endl; std::cout << "pen color_surface = lightred;" <<
  //   std::endl; std::cout << "material material_face =
  //   material(diffusepen=gray(0.8),
  //   "
  //                "emissivepen=gray(0.6), specularpen=black);"
  //             << std::endl;
  //   std::cout << "material material_surface =
  //   material(lightred);" << std::endl; std::cout << "bool
  //   plot_cell = true;" << std::endl; std::cout
  //   << "bool plot_face_label = true;" << std::endl; std::cout
  //   << "plot_face_normals = true;" << std::endl; std::cout <<
  //   "bool plot_cell_unclipped = true;" << std::endl; std::cout
  //   << "bool plot_face_unclipped_label = true;" << std::endl;
  //   std::cout << "bool plot_proj_cell_unclipped = true;" <<
  //   std::endl; std::cout << "bool plot_cell_edges = true;" <<
  //   std::endl; std::cout << "bool plot_cell_unclipped_edges =
  //   true;" << std::endl; std::cout << "bool
  //   plot_proj_cell_unclipped_edges = true;" << std::endl;
  //   std::cout << "bool plot_clipped_corners = true;" <<
  //   std::endl; std::cout << "bool plot_unclipped_corners =
  //   true;" << std::endl; std::cout << "bool plot_intersections
  //   = true;" << std::endl; std::cout << "bool plot_arcs =
  //   true;" << std::endl; std::cout << "bool plot_proj_arcs =
  //   true;" << std::endl; std::cout << "bool plot_surface =
  //   true;" << std::endl; std::cout << "bool plot_paraboloid =
  //   true;" << std::endl; std::cout << "bool
  //   plot_paraboloid_label = true;" << std::endl; std::cout <<
  //   "real link_arrow_size = 0;" << std::endl; std::cout <<
  //   "real link_size = 0.25;" << std::endl; std::cout << "real
  //   proj_x = -0.6;" << std::endl; std::cout << "real proj_y =
  //   -0.8;" << std::endl; std::cout << "real proj_z = -1.25;" <<
  //   std::endl; std::cout << "real R = 1.0;" << std::endl;
  //   std::cout << "real normal_length = 0.5;" << std::endl;

  //   ////////////////// VERTICES
  //   std::cout << std::endl
  //             << "// " << vertices_after_intersection << "
  //             vertices"
  //             << std::endl;
  //   for (UnsignedIndex_t v = 0; v <
  //   vertices_after_intersection; ++v) {
  //     const Pt& pt =
  //     a_polytope->getVertex(v)->getLocation().getPt();
  //     std::cout << std::setprecision(20) << "triple P"
  //               << reinterpret_cast<std::uintptr_t>(&pt) << " =
  //               (" << pt[0]
  //               << "," << pt[1] << "," << pt[2] << ");" <<
  //               std::endl;
  //     std::cout << std::setprecision(20) << "triple Pxy"
  //               << reinterpret_cast<std::uintptr_t>(&pt) << " =
  //               (" << pt[0]
  //               << "," << pt[1] << ",proj_z);" << std::endl;
  //   }

  //   ////////////////////////// FACES
  //   std::cout << "// " << starting_number_of_faces << " faces"
  //   << std::endl; for (UnsignedIndex_t f = 0; f <
  //   starting_number_of_faces;
  //   ++f) {
  //     auto center = Pt(0.0, 0.0, 0.0);
  //     auto tmp_center = Pt(0.0, 0.0, 0.0);
  //     auto& face = *(*a_polytope)[f];
  //     const auto starting_half_edge =
  //     face.getStartingHalfEdge(); auto current_half_edge =
  //     starting_half_edge; std::cout << "path3 face" << f <<
  //     "_polygon = "; int count_vert = 0; do {
  //       const auto& vertex = *(current_half_edge->getVertex());
  //       const Pt& pt = vertex.getLocation().getPt();
  //       std::cout << "P" <<
  //       reinterpret_cast<std::uintptr_t>(&pt) << " --
  //       "; current_half_edge =
  //       current_half_edge->getNextHalfEdge(); if
  //       (vertex.needsToSeek()) {
  //         center += pt;
  //         count_vert++;
  //         if (count_vert >= 2 && count_vert <= 4) {
  //           tmp_center += pt;
  //         }
  //       }
  //     } while (current_half_edge != starting_half_edge);
  //     center /= static_cast<double>(count_vert);
  //     std::cout << "cycle;" << std::endl;
  //     std::cout << "surface face" << f << " = surface(face" <<
  //     f <<
  //     "_polygon);"
  //               << std::endl;
  //     std::cout << std::setprecision(20) << "triple C" << f <<
  //     " = ("
  //               << center[0] << "," << center[1] << "," <<
  //               center[2] <<
  //               ");"
  //               << std::endl;
  //     std::cout << "if (plot_face_label) "
  //                  "{label(scale3(1.25)*align(rotate(00)*\"$\\mathcal{F}_{"
  //               << f + 1 << "}$\",normal(face" << f <<
  //               "_polygon)),C" <<
  //               f
  //               << "+0.001*normal(face" << f <<
  //               "_polygon),black,Embedded);}"
  //               << std::endl;
  //     tmp_center /= 3.0;
  //     std::cout << std::setprecision(20) << "triple TmpC" << f
  //     << " = ("
  //               << tmp_center[0] << "," << tmp_center[1] << ","
  //               << tmp_center[2]
  //               << ");" << std::endl;
  //     std::cout << "triple tmp_end_point" << f << " = TmpC" <<
  //     f
  //               << " + normal_length * normal(face" << f <<
  //               "_polygon);"
  //               << std::endl;
  //     std::cout << "if (plot_face_normals) "
  //                  "{path3 patharrowF"
  //               << f << " =TmpC" << f << " -- tmp_end_point" <<
  //               f
  //               << "; "
  //                  "label(scale3(1.25)*align(rotate(00)*\"\\tiny$\\times$\","
  //                  "normal(face"
  //               << f << "_polygon)),TmpC" << f <<
  //               "+0.001*normal(face" <<
  //               f
  //               << "_polygon),black,Embedded); draw(patharrowF"
  //               << f
  //               << ", arrow = Arrow3(emissive(black), size =
  //               4bp),black +
  //               "
  //                  "linewidth(0.5pt),L =
  //                  Label(\"\\scriptsize$\\mathbf{n}_{"
  //               << f + 1 << "}$\", position = EndPoint, align =
  //               E));}"
  //               << std::endl;
  //   }
  //   std::cout << "surface poly = surface(face";
  //   for (UnsignedIndex_t f = 0; f < starting_number_of_faces -
  //   1; ++f) {
  //     std::cout << f << ",face";
  //   }
  //   std::cout << starting_number_of_faces - 1 << ");" <<
  //   std::endl;

  //   ////////////////////////// ARCS AND SURFACE
  //   const std::vector<RationalBezierArc>& arcs =
  //   a_surface->getArcs(); UnsignedIndex_t nArcs = arcs.size();
  //   std::cout << "// " << nArcs << " arcs" << std::endl;
  //   Pt center_surf = Pt(0.0, 0.0, 0.0);
  //   for (UnsignedIndex_t a = 0; a < nArcs; ++a) {
  //     std::cout << std::setprecision(20) << "real W" << a << "
  //     = "
  //               << arcs[a].weight() << ";" << std::endl;
  //     center_surf += arcs[a].start_point().getPt();
  //   }
  //   center_surf /= static_cast<double>(nArcs);
  //   center_surf[2] =
  //       -a_aligned_paraboloid.a() * center_surf[0] *
  //       center_surf[0] - a_aligned_paraboloid.b() *
  //       center_surf[1] * center_surf[1];
  //   std::cout << std::setprecision(20) << "triple P0 = (" <<
  //   center_surf[0]
  //             << "," << center_surf[1] << "," << center_surf[2]
  //             << ");"
  //             << std::endl;
  //   for (UnsignedIndex_t a = 0; a < nArcs; ++a) {
  //     const Pt& start_pt = arcs[a].start_point();
  //     const Pt& control_pt = arcs[a].control_point();
  //     const Pt& end_pt = arcs[a].end_point();
  //     std::cout << std::setprecision(20) << "triple CP" << a <<
  //     " = ("
  //               << control_pt[0] << "," << control_pt[1] << ","
  //               << control_pt[2]
  //               << ");" << std::endl;
  //     std::cout << std::setprecision(20) << "triple CPxy" << a
  //     << " = ("
  //               << control_pt[0] << "," << control_pt[1] <<
  //               ",proj_z);"
  //               << std::endl;
  //     std::cout << "triple arc"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "f(real t) { return ((P"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "*(1-t)^2 + 2*W" << a << "*CP" << a <<
  //               "*(1-t)*t + P"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "*t^2) / ((1-t)^2 + 2*W" << a << "*(1-t)*t +
  //               t^2));}"
  //               << std::endl;
  //     std::cout << "triple arcxy"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "f(real t) { return ((Pxy"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "*(1-t)^2 + 2*W" << a << "*CPxy" << a <<
  //               "*(1-t)*t + Pxy"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "*t^2) / ((1-t)^2 + 2*W" << a << "*(1-t)*t +
  //               t^2));}"
  //               << std::endl;
  //     std::cout << "path3 arc"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << " = graph(arc"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "f, 0, 1, operator ..);" << std::endl;
  //     std::cout << "path3 arcxy"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << " = graph(arcxy"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "f, 0, 1, operator ..);" << std::endl;
  //     std::cout << "triple arc"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "f(real t) { return ((P"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "*(1-t)^2 + 2*W" << a << "*CP" << a <<
  //               "*(1-t)*t + P"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "*t^2) / ((1-t)^2 + 2*W" << a << "*(1-t)*t +
  //               t^2));}"
  //               << std::endl;
  //     std::cout << "triple arcxy"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "f(real t) { return ((Pxy"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "*(1-t)^2 + 2*W" << a << "*CPxy" << a <<
  //               "*(1-t)*t + Pxy"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "*t^2) / ((1-t)^2 + 2*W" << a << "*(1-t)*t +
  //               t^2));}"
  //               << std::endl;
  //     std::cout << "path3 arc"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << " = graph(arc"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "f, 0, 1, operator ..);" << std::endl;
  //     std::cout << "path3 arcxy"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << " = graph(arcxy"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "f, 0, 1, operator ..);" << std::endl;
  //     std::cout << "triple P"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "P0(real t) {triple T = P0*t+P"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "*(1-t); return(T.x,T.y,(" <<
  //               std::setprecision(20)
  //               << -a_aligned_paraboloid.a() << ")*T.x^2+("
  //               << -a_aligned_paraboloid.b() << ")*T.y^2);}" <<
  //               std::endl;
  //   }
  //   for (UnsignedIndex_t a = 0; a < nArcs; ++a) {
  //     std::cout << "triple tri0to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "_func(pair t) { real u = t.x, v = t.y;"
  //                  "triple T = arc"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "f(u) * (1 - v) + P0 * v + "
  //                  "P"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "P0(v) * (1 - u) + P"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "* u - P"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << " * (1 - u) * (1 - v) - P0 * (1 - u) * v -
  //               P"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << " * u * (1 - v) - P0 * u * v;"
  //                  "return (T.x, T.y, ("
  //               << -a_aligned_paraboloid.a() << ")*T.x^2 +("
  //               << -a_aligned_paraboloid.b() << ")*T.y^2);}" <<
  //               std::endl;
  //     std::cout << "surface tri0to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << " = surface(tri0to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << "_func,(0,0),(1,1),100,100,Spline);" <<
  //               std::endl;
  //     std::cout << "if (plot_surface) {draw(tri0to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << ", material_surface);} " << std::endl;
  //   }

  //   ////////////////////////// UNCLIPPED FACES
  //   std::cout << "// " << starting_number_of_faces << " faces"
  //   << std::endl; int
  //   unclipped_faces[starting_number_of_faces]; for
  //   (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f)
  //   {
  //     auto& face = *(*a_polytope)[f];
  //     const auto intersection_size =
  //     face.getNumberOfIntersections(); const auto
  //     starting_half_edge = face.getStartingHalfEdge(); if
  //     (intersection_size > 0 ||
  //         starting_half_edge->getVertex()->isNotClipped()) {
  //       unclipped_faces[f] = 1;
  //       auto current_half_edge = starting_half_edge;
  //       std::cout << "path3 face" << f << "_unclipped_polygon =
  //       "; bool face_exists = false; do {
  //         const auto& vertex =
  //         *(current_half_edge->getVertex()); const auto&
  //         next_vertex =
  //             *(current_half_edge->getNextHalfEdge()->getVertex());
  //         if (vertex.needsToSeek()) {
  //           const Pt& pt = vertex.getLocation().getPt();
  //           double distance = signedDistance(pt,
  //           a_aligned_paraboloid); if (distance < -TEN *
  //           DBL_EPSILON) {
  //             std::cout << "P" <<
  //             reinterpret_cast<std::uintptr_t>(&pt)
  //                       << " -- ";
  //             face_exists = true;
  //             if (next_vertex.doesNotNeedToSeek())  //
  //             intersection
  //             {
  //               current_half_edge =
  //               current_half_edge->getNextHalfEdge(); const Pt&
  //               next_pt = next_vertex.getLocation().getPt();
  //               std::cout << "P" <<
  //               reinterpret_cast<std::uintptr_t>(&next_pt)
  //                         << " -- ";
  //               // Find arc
  //               std::uintptr_t candidates[2];
  //               int arc_found = 0;
  //               for (UnsignedIndex_t a = 0; a < nArcs; ++a) {
  //                 const Pt& start_pt = arcs[a].start_point();
  //                 const Pt& control_pt =
  //                 arcs[a].control_point(); const Pt& end_pt =
  //                 arcs[a].end_point(); if
  //                 (arcs[a].start_point_id() ==
  //                     reinterpret_cast<std::uintptr_t>(&next_pt))
  //                     {
  //                   candidates[0] = arcs[a].end_point_id();
  //                   arc_found++;
  //                 } else if (arcs[a].end_point_id() ==
  //                            reinterpret_cast<std::uintptr_t>(&next_pt))
  //                            {
  //                   candidates[1] = arcs[a].start_point_id();
  //                   arc_found++;
  //                 }
  //               }
  //               if (arc_found != 2) exit(-1);
  //               bool next_intersection_found = false;
  //               while (!next_intersection_found) {
  //                 current_half_edge =
  //                 current_half_edge->getNextHalfEdge(); const
  //                 auto& test_vertex =
  //                 *(current_half_edge->getVertex()); const Pt&
  //                 test_pt = test_vertex.getLocation().getPt();
  //                 if (candidates[0] ==
  //                     reinterpret_cast<std::uintptr_t>(&test_pt))
  //                     {
  //                   next_intersection_found = true;
  //                 } else if (candidates[1] ==
  //                            reinterpret_cast<std::uintptr_t>(&test_pt))
  //                            {
  //                   next_intersection_found = true;
  //                 }
  //                 if (next_intersection_found) {
  //                   std::cout
  //                       << "arc" <<
  //                       reinterpret_cast<std::uintptr_t>(&next_pt)
  //                       << "to" <<
  //                       reinterpret_cast<std::uintptr_t>(&test_pt)
  //                       << " -- ";
  //                 }
  //               }
  //               if (current_half_edge == starting_half_edge) {
  //                 break;
  //               }
  //             }
  //           }
  //         }
  //         current_half_edge =
  //         current_half_edge->getNextHalfEdge();
  //       } while (current_half_edge != starting_half_edge);
  //       if (face_exists) {
  //         std::cout << "cycle;" << std::endl;
  //       } else {
  //         std::cout << ";" << std::endl;
  //       }
  //       std::cout << "surface face_unclipped" << f << " =
  //       surface(face"
  //       <<
  //       f
  //                 << "_unclipped_polygon);" << std::endl;
  //       std::cout
  //           << "if (plot_face_unclipped_label) "
  //              "{label(scale3(1.25)*align(rotate(00)*\"$\\hat{\\mathcal{F}}_{"
  //           << f + 1 << "}$\",normal(face" << f <<
  //           "_polygon)),C" << f
  //           << "+0.001*normal(face" << f <<
  //           "_polygon),black,Embedded);}"
  //           << std::endl;

  //     } else {
  //       unclipped_faces[f] = 0;
  //     }
  //   }
  //   std::cout << "surface poly_unclipped = surface(";
  //   bool comma = false;
  //   for (UnsignedIndex_t f = 0; f < starting_number_of_faces;
  //   ++f) {
  //     if (unclipped_faces[f]) {
  //       if (comma) {
  //         std::cout << ",face_unclipped" << f;
  //       } else {
  //         std::cout << "face_unclipped" << f;
  //         comma = true;
  //       }
  //     }
  //   }
  //   std::cout << ");" << std::endl;

  //   ///////////////////////////////// PROJ UNCLIPPED CELL
  //   for (UnsignedIndex_t f = 0; f < starting_number_of_faces;
  //   ++f) {
  //     auto& face = *(*a_polytope)[f];
  //     const auto intersection_size =
  //     face.getNumberOfIntersections(); const auto
  //     starting_half_edge = face.getStartingHalfEdge(); if
  //     (intersection_size > 0 ||
  //         starting_half_edge->getVertex()->isNotClipped()) {
  //       unclipped_faces[f] = 1;
  //       auto current_half_edge = starting_half_edge;
  //       std::cout << "path3 face" << f <<
  //       "_proj_unclipped_polygon = "; bool face_exists = false;
  //       do {
  //         const auto& vertex =
  //         *(current_half_edge->getVertex()); const auto&
  //         next_vertex =
  //             *(current_half_edge->getNextHalfEdge()->getVertex());
  //         if (vertex.needsToSeek()) {
  //           const Pt& pt = vertex.getLocation().getPt();
  //           double distance = signedDistance(pt,
  //           a_aligned_paraboloid); if (distance < -TEN *
  //           DBL_EPSILON) {
  //             std::cout << "Pxy" <<
  //             reinterpret_cast<std::uintptr_t>(&pt)
  //                       << " -- ";
  //             face_exists = true;
  //             if (next_vertex.doesNotNeedToSeek())  //
  //             intersection
  //             {
  //               current_half_edge =
  //               current_half_edge->getNextHalfEdge(); const Pt&
  //               next_pt = next_vertex.getLocation().getPt();
  //               std::cout << "Pxy" <<
  //               reinterpret_cast<std::uintptr_t>(&next_pt)
  //                         << " -- ";
  //               // Find arc
  //               std::uintptr_t candidates[2];
  //               int arc_found = 0;
  //               for (UnsignedIndex_t a = 0; a < nArcs; ++a) {
  //                 const Pt& start_pt = arcs[a].start_point();
  //                 const Pt& control_pt =
  //                 arcs[a].control_point(); const Pt& end_pt =
  //                 arcs[a].end_point(); if
  //                 (arcs[a].start_point_id() ==
  //                     reinterpret_cast<std::uintptr_t>(&next_pt))
  //                     {
  //                   candidates[0] = arcs[a].end_point_id();
  //                   arc_found++;
  //                 } else if (arcs[a].end_point_id() ==
  //                            reinterpret_cast<std::uintptr_t>(&next_pt))
  //                            {
  //                   candidates[1] = arcs[a].start_point_id();
  //                   arc_found++;
  //                 }
  //               }
  //               if (arc_found != 2) exit(-1);
  //               bool next_intersection_found = false;
  //               while (!next_intersection_found) {
  //                 current_half_edge =
  //                 current_half_edge->getNextHalfEdge(); const
  //                 auto& test_vertex =
  //                 *(current_half_edge->getVertex()); const Pt&
  //                 test_pt = test_vertex.getLocation().getPt();
  //                 if (candidates[0] ==
  //                     reinterpret_cast<std::uintptr_t>(&test_pt))
  //                     {
  //                   next_intersection_found = true;
  //                 } else if (candidates[1] ==
  //                            reinterpret_cast<std::uintptr_t>(&test_pt))
  //                            {
  //                   next_intersection_found = true;
  //                 }
  //                 if (next_intersection_found) {
  //                   std::cout
  //                       << "arcxy" <<
  //                       reinterpret_cast<std::uintptr_t>(&next_pt)
  //                       << "to" <<
  //                       reinterpret_cast<std::uintptr_t>(&test_pt)
  //                       << " -- ";
  //                 }
  //               }
  //               if (current_half_edge == starting_half_edge) {
  //                 break;
  //               }
  //             }
  //           }
  //         }
  //         current_half_edge =
  //         current_half_edge->getNextHalfEdge();
  //       } while (current_half_edge != starting_half_edge);
  //       if (face_exists) {
  //         std::cout << "cycle;" << std::endl;
  //       } else {
  //         std::cout << ";" << std::endl;
  //       }
  //       std::cout << "surface face_proj_unclipped" << f << " =
  //       surface(face"
  //                 << f << "_proj_unclipped_polygon);" <<
  //                 std::endl;
  //     } else {
  //       unclipped_faces[f] = 0;
  //     }
  //   }
  //   std::cout << "surface poly_proj_unclipped = surface(";
  //   comma = false;
  //   for (UnsignedIndex_t f = 0; f < starting_number_of_faces;
  //   ++f) {
  //     if (unclipped_faces[f]) {
  //       if (comma) {
  //         std::cout << ",face_proj_unclipped" << f;
  //       } else {
  //         std::cout << "face_proj_unclipped" << f;
  //         comma = true;
  //       }
  //     }
  //   }
  //   std::cout << ");" << std::endl;

  //   for (UnsignedIndex_t v = 0; v <
  //   vertices_after_intersection; ++v) {
  //     const Pt& pt =
  //     a_polytope->getVertex(v)->getLocation().getPt(); double
  //     distance = signedDistance(pt, a_aligned_paraboloid); if
  //     (distance < -TEN * DBL_EPSILON) {
  //       std::cout
  //           << "if (plot_unclipped_corners) {draw(shift(P"
  //           << reinterpret_cast<std::uintptr_t>(&pt)
  //           << ")*scale3(scale_vertex)*unitsphere,
  //           color_corner_unclipped);}"
  //           << std::endl;
  //       std::cout << "if (plot_proj_links) {draw(P"
  //                 << reinterpret_cast<std::uintptr_t>(&pt) <<
  //                 "--Pxy"
  //                 << reinterpret_cast<std::uintptr_t>(&pt)
  //                 << ", "
  //                    "arrow=MidArrow3(emissive(color_unclipped_edge),size=link_"
  //                    "arrow_size), "
  //                    "color_unclipped_edge+linewidth(link_size));}"
  //                 << std::endl;
  //     } else if (distance > TEN * DBL_EPSILON) {
  //       std::cout
  //           << "if (plot_clipped_corners) {draw(shift(P"
  //           << reinterpret_cast<std::uintptr_t>(&pt)
  //           << ")*scale3(scale_vertex)*unitsphere,
  //           color_corner_clipped);}"
  //           << std::endl;
  //     } else {
  //       std::cout << "if (plot_intersections) {draw(shift(P"
  //                 << reinterpret_cast<std::uintptr_t>(&pt)
  //                 << ")*scale3(scale_vertex)*unitsphere,
  //                 color_intersection);}"
  //                 << std::endl;
  //       std::cout
  //           << "if (plot_proj_links) {draw(P"
  //           << reinterpret_cast<std::uintptr_t>(&pt) << "--Pxy"
  //           << reinterpret_cast<std::uintptr_t>(&pt)
  //           << ",
  //           arrow=MidArrow3(emissive(color_arc),size=link_arrow_size),
  //           "
  //              "color_arc+linewidth(link_size));}"
  //           << std::endl;
  //     }
  //   }
  //   std::cout << "transform T_tube(real t){return
  //   scale(scale_edge);}"
  //             << std::endl;
  //   std::cout << "transform T_tube_unclipped(real t){return "
  //                "scale(scale_unclipped_edge);}"
  //             << std::endl;
  //   std::cout << "transform T_arc(real t){return
  //   scale(scale_arc);}"
  //             << std::endl;
  //   std::cout << "if (plot_cell_edges) {" << std::endl;
  //   for (UnsignedIndex_t f = 0; f < starting_number_of_faces -
  //   1; ++f) {
  //     std::cout << "draw(tube(face" << f
  //               << "_polygon, unitcircle, T_tube),
  //               color_edge);" << std::endl;
  //   }
  //   std::cout << "}" << std::endl;
  //   std::cout << "if (plot_cell_unclipped_edges) {" <<
  //   std::endl; for (UnsignedIndex_t f = 0; f <
  //   starting_number_of_faces - 1; ++f) {
  //     if (unclipped_faces[f]) {
  //       std::cout << "draw(tube(face" << f
  //                 << "_unclipped_polygon, unitcircle,
  //                 T_tube_unclipped),
  //                 "
  //                    "color_unclipped_edge);"
  //                 << std::endl;
  //     }
  //   }
  //   std::cout << "}" << std::endl;
  //   std::cout << "if (plot_proj_cell_unclipped_edges) {" <<
  //   std::endl; for (UnsignedIndex_t f = 0; f <
  //   starting_number_of_faces - 1; ++f) {
  //     if (unclipped_faces[f]) {
  //       std::cout << "draw(face" << f
  //                 << "_proj_unclipped_polygon, "
  //                    "color_unclipped_edge+linewidth(0.75pt));"
  //                 << std::endl;
  //     }
  //   }
  //   std::cout << "}" << std::endl;

  //   std::cout << "if (plot_arcs) {" << std::endl;
  //   for (UnsignedIndex_t a = 0; a < nArcs; ++a) {
  //     std::cout << "draw(tube(arc"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << ", unitcircle, T_arc), color_arc);" <<
  //               std::endl;
  //   }
  //   std::cout << "}" << std::endl;
  //   std::cout << "if (plot_proj_arcs) {" << std::endl;
  //   for (UnsignedIndex_t a = 0; a < nArcs; ++a) {
  //     std::cout << "draw(arcxy"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].start_point_id())
  //               << "to"
  //               <<
  //               reinterpret_cast<std::uintptr_t>(arcs[a].end_point_id())
  //               << ", color_arc+linewidth(1pt));" << std::endl;
  //   }
  //   std::cout << "}" << std::endl;
  //   std::cout << "if (plot_cell) {draw(poly, material_face);}"
  //   << std::endl; std::cout
  //       << "if (plot_cell_unclipped) {draw(poly_unclipped,
  //       material_face);}"
  //       << std::endl;
  //   std::cout << "if (plot_surface_label) {real gamma_x = 0.2,
  //   gamma_y = 0.1;
  //   "
  //                "triple gamma_loc = (gamma_x,gamma_y,"
  //             << -a_aligned_paraboloid.a() << "*gamma_x^2 + "
  //             << -a_aligned_paraboloid.b()
  //             << "*gamma_y^2); triple normal_gamma_loc = (2.0 *
  //             "
  //             << a_aligned_paraboloid.a() << "*gamma_x, 2.0 * "
  //             << a_aligned_paraboloid.b()
  //             << "*gamma_y, 1); label(scale3(1.25) *
  //             align(rotate(-20) *
  //             "
  //                "\"$\\tilde{\\mathcal{S}}$\",
  //                normal_gamma_loc), gamma_loc
  //                + " "0.01 * normal_gamma_loc, color_arc,
  //                Embedded);}"
  //             << std::endl;
  //   std::cout << "if (plot_paraboloid) {triple
  //   paraboloid_func(pair t) { real
  //   "
  //                "u = t.x*cos(t.y), v = t.x*sin(t.y);"
  //                "return (u, v, ("
  //             << -a_aligned_paraboloid.a() << ")*u^2 +("
  //             << -a_aligned_paraboloid.b() << ")*v^2);}" <<
  //             std::endl;
  //   std::cout << "surface paraboloid = "
  //                "surface(paraboloid_func,(0,0),(R,2*pi),200,200,Spline);"
  //             << std::endl;
  //   std::cout << "draw(paraboloid, material_surface);} " <<
  //   std::endl; std::cout
  //       << "if (plot_paraboloid_label) {real gamma_x = 1.0,
  //       gamma_y = 0.0;
  //       "
  //          "triple gamma_loc = (gamma_x,gamma_y,"
  //       << -a_aligned_paraboloid.a() << "*gamma_x^2 + "
  //       << -a_aligned_paraboloid.b()
  //       << "*gamma_y^2); triple normal_gamma_loc = (2.0 * "
  //       << a_aligned_paraboloid.a() << "*gamma_x, 2.0 * "
  //       << a_aligned_paraboloid.b()
  //       << "*gamma_y, 1); label(scale3(1.25) *
  //       align(rotate(-20) * "
  //          "\"${\\mathcal{S}}$\", normal_gamma_loc), gamma_loc
  //          + " "0.01 * normal_gamma_loc, color_arc,
  //          Embedded);}"
  //       << std::endl;
  //   std::cout << "/////// END ASYMPTOTE CODE" << std::endl
  //             << std::endl
  //             << std::endl;

  //   std::cout << "%%%%%%%%%% START TIKZ CODE" << std::endl;
  //   char tikz_colors[12][1000] = {"Red3",        "Green3",
  //   "Blue4",
  //                                 "DarkOrchid3", "Coral3",
  //                                 "LightBlue4", "Tomato4",
  //                                 "DeepPink2", "OliveDrab4",
  //                                 "Azure4", "DeepSkyBlue2",
  //                                 "DarkOrchid1"};
  //   ////////////////// VERTICES
  //   std::cout << std::endl
  //             << "%% " << vertices_after_intersection << "
  //             vertices"
  //             << std::endl;
  //   for (UnsignedIndex_t v = 0; v <
  //   vertices_after_intersection; ++v) {
  //     const Pt& pt =
  //     a_polytope->getVertex(v)->getLocation().getPt();
  //     std::cout << std::setprecision(20) << "\\coordinate (Pxy"
  //               << reinterpret_cast<std::uintptr_t>(&pt) << ")
  //               at (" << pt[0]
  //               << "," << pt[1] << ");" << std::endl;
  //   }
  //   ////////////////////////// ARCS AND SURFACE
  //   std::cout << "%% " << nArcs << " arcs" << std::endl;
  //   char arcs_string[2 * nArcs][1000];
  //   char arcs_arrow_string[2 * nArcs][1000];
  //   for (UnsignedIndex_t a = 0; a < nArcs; ++a) {
  //     const Pt& start_pt = arcs[a].start_point();
  //     const Pt& control_pt = arcs[a].control_point();
  //     const Pt& end_pt = arcs[a].end_point();
  //     // std::cout << std::setprecision(20) << "\\coordinate
  //     (CPxy" << a
  //     //           << ") at (" << control_pt[0] << "," <<
  //     control_pt[1]
  //     <<
  //     ");"
  //     //           << std::endl;
  //     sprintf(arcs_string[a],
  //             "plot[domain=0:1,variable=\\x,samples=100] "
  //             "({(%.15e*(1.0-\\x)*(1.0-\\x)
  //             + 2.0*%.15e*%.15e*(1.0-\\x)*\\x
  //             +
  //             "
  //             "%.15e*\\x*\\x) / ((1.0-\\x)*(1.0-\\x)
  //             + 2.0*%.15e*(1.0-\\x)*\\x "
  //             "+ \\x*\\x)},{(%.15e*(1.0-\\x)*(1.0-\\x) + "
  //             "2.0*%.15e*%.15e*(1.0-\\x)*\\x + %.15e*\\x*\\x) /
  //             "
  //             "((1.0-\\x)*(1.0-\\x) + 2.0*%.15e*(1.0-\\x)*\\x +
  //             \\x*\\x)})", start_pt[0], arcs[a].weight(),
  //             control_pt[0], end_pt[0], arcs[a].weight(),
  //             start_pt[1], arcs[a].weight(), control_pt[1],
  //             end_pt[1], arcs[a].weight());
  //     sprintf(arcs_arrow_string[a],
  //             "plot[domain=0:0.5,variable=\\x,samples=100] "
  //             "({(%.15e*(1.0-\\x)*(1.0-\\x)
  //             + 2.0*%.15e*%.15e*(1.0-\\x)*\\x
  //             +
  //             "
  //             "%.15e*\\x*\\x) / ((1.0-\\x)*(1.0-\\x)
  //             + 2.0*%.15e*(1.0-\\x)*\\x "
  //             "+ \\x*\\x)},{(%.15e*(1.0-\\x)*(1.0-\\x) + "
  //             "2.0*%.15e*%.15e*(1.0-\\x)*\\x + %.15e*\\x*\\x) /
  //             "
  //             "((1.0-\\x)*(1.0-\\x) + 2.0*%.15e*(1.0-\\x)*\\x +
  //             \\x*\\x)})", start_pt[0], arcs[a].weight(),
  //             control_pt[0], end_pt[0], arcs[a].weight(),
  //             start_pt[1], arcs[a].weight(), control_pt[1],
  //             end_pt[1], arcs[a].weight());
  //     sprintf(arcs_string[nArcs + a],
  //             "plot[domain=0:1,variable=\\x,samples=100] "
  //             "({(%.15e*(1.0-\\x)*(1.0-\\x)
  //             + 2.0*%.15e*%.15e*(1.0-\\x)*\\x
  //             +
  //             "
  //             "%.15e*\\x*\\x) / ((1.0-\\x)*(1.0-\\x)
  //             + 2.0*%.15e*(1.0-\\x)*\\x "
  //             "+ \\x*\\x)},{(%.15e*(1.0-\\x)*(1.0-\\x) + "
  //             "2.0*%.15e*%.15e*(1.0-\\x)*\\x + %.15e*\\x*\\x) /
  //             "
  //             "((1.0-\\x)*(1.0-\\x) + 2.0*%.15e*(1.0-\\x)*\\x +
  //             \\x*\\x)})", end_pt[0], arcs[a].weight(),
  //             control_pt[0], start_pt[0], arcs[a].weight(),
  //             end_pt[1], arcs[a].weight(), control_pt[1],
  //             start_pt[1], arcs[a].weight());
  //     sprintf(arcs_arrow_string[nArcs + a],
  //             "plot[domain=0:0.5,variable=\\x,samples=100] "
  //             "({(%.15e*(1.0-\\x)*(1.0-\\x)
  //             + 2.0*%.15e*%.15e*(1.0-\\x)*\\x
  //             +
  //             "
  //             "%.15e*\\x*\\x) / ((1.0-\\x)*(1.0-\\x)
  //             + 2.0*%.15e*(1.0-\\x)*\\x "
  //             "+ \\x*\\x)},{(%.15e*(1.0-\\x)*(1.0-\\x) + "
  //             "2.0*%.15e*%.15e*(1.0-\\x)*\\x + %.15e*\\x*\\x) /
  //             "
  //             "((1.0-\\x)*(1.0-\\x) + 2.0*%.15e*(1.0-\\x)*\\x +
  //             \\x*\\x)})", end_pt[0], arcs[a].weight(),
  //             control_pt[0], start_pt[0], arcs[a].weight(),
  //             end_pt[1], arcs[a].weight(), control_pt[1],
  //             start_pt[1], arcs[a].weight());
  //   }
  //   ///////////////////////////////// PROJ UNCLIPPED CELL
  //   char faces_string[starting_number_of_faces][1000];
  //   char faces_arrow_string[starting_number_of_faces][1000];
  //   for (UnsignedIndex_t f = 0; f < starting_number_of_faces;
  //   ++f) {
  //     auto& face = *(*a_polytope)[f];
  //     auto face_normal = face.getPlane().normal();
  //     const auto intersection_size =
  //     face.getNumberOfIntersections(); const auto
  //     starting_half_edge = face.getStartingHalfEdge(); auto
  //     center = Pt(0.0, 0.0, 0.0); double scale_center = 0.0; if
  //     (intersection_size > 0 ||
  //         starting_half_edge->getVertex()->isNotClipped()) {
  //       std::cout << "%% FACE " << f << std::endl;
  //       std::cout << (face_normal[2] > 0.0 ?
  //       "\\begin{pgfonlayer}{foreground}"
  //                                          :
  //                                          "\\begin{pgfonlayer}{background}")
  //                 << std::endl;
  //       std::cout << "\\begin{scope}" << std::endl;
  //       unclipped_faces[f] = 1;
  //       auto current_half_edge = starting_half_edge;
  //       sprintf(faces_string[f], " ");
  //       sprintf(faces_arrow_string[f], " ");
  //       bool face_exists = false;
  //       do {
  //         const auto& vertex =
  //         *(current_half_edge->getVertex()); const auto&
  //         next_vertex =
  //             *(current_half_edge->getNextHalfEdge()->getVertex());
  //         if (vertex.needsToSeek()) {
  //           const Pt& pt = vertex.getLocation().getPt();
  //           double distance = signedDistance(pt,
  //           a_aligned_paraboloid); if (distance < -TEN *
  //           DBL_EPSILON) {
  //             sprintf(faces_string[f], "%s(Pxy%ld) -- ",
  //             faces_string[f],
  //                     reinterpret_cast<std::uintptr_t>(&pt));
  //             center += pt;
  //             scale_center += 1.0;
  //             sprintf(faces_arrow_string[f],
  //                     "%s\\draw[very thick,-{angle
  //                     60[scale=0.75]},%s] "
  //                     "(Pxy%ld) -- "
  //                     "(%.15e,%.15e);\n",
  //                     faces_arrow_string[f],
  //                     face_normal[2] > 0.0 ? "Red3" : "Blue4",
  //                     reinterpret_cast<std::uintptr_t>(&pt),
  //                     0.5 * (pt[0] +
  //                     next_vertex.getLocation().getPt()[0]),
  //                     0.5 * (pt[1]
  //                     + next_vertex.getLocation().getPt()[1]));
  //             face_exists = true;
  //             if (next_vertex.doesNotNeedToSeek())  //
  //             intersection
  //             {
  //               current_half_edge =
  //               current_half_edge->getNextHalfEdge(); const Pt&
  //               next_pt = next_vertex.getLocation().getPt();
  //               // sprintf(faces_string[f], "%s(Pxy%ld) -- ",
  //               faces_string[f],
  //               // reinterpret_cast<std::uintptr_t>(&next_pt));
  //               // Find arc
  //               center += next_pt;
  //               scale_center += 1.0;
  //               std::uintptr_t candidates[2];
  //               int arc_found = 0;
  //               for (UnsignedIndex_t a = 0; a < nArcs; ++a) {
  //                 const Pt& start_pt = arcs[a].start_point();
  //                 const Pt& control_pt =
  //                 arcs[a].control_point(); const Pt& end_pt =
  //                 arcs[a].end_point(); if
  //                 (arcs[a].start_point_id() ==
  //                     reinterpret_cast<std::uintptr_t>(&next_pt))
  //                     {
  //                   candidates[0] = arcs[a].end_point_id();
  //                   arc_found++;
  //                 } else if (arcs[a].end_point_id() ==
  //                            reinterpret_cast<std::uintptr_t>(&next_pt))
  //                            {
  //                   candidates[1] = arcs[a].start_point_id();
  //                   arc_found++;
  //                 }
  //               }
  //               if (arc_found != 2) exit(-1);
  //               bool next_intersection_found = false;
  //               while (!next_intersection_found) {
  //                 current_half_edge =
  //                 current_half_edge->getNextHalfEdge(); const
  //                 auto& test_vertex =
  //                 *(current_half_edge->getVertex()); const Pt&
  //                 test_pt = test_vertex.getLocation().getPt();
  //                 if (candidates[0] ==
  //                     reinterpret_cast<std::uintptr_t>(&test_pt))
  //                     {
  //                   next_intersection_found = true;
  //                 } else if (candidates[1] ==
  //                            reinterpret_cast<std::uintptr_t>(&test_pt))
  //                            {
  //                   next_intersection_found = true;
  //                 }
  //                 if (next_intersection_found) {
  //                   for (UnsignedIndex_t a = 0; a < nArcs; ++a)
  //                   {
  //                     if (arcs[a].start_point_id() ==
  //                             reinterpret_cast<std::uintptr_t>(&next_pt)
  //                             &&
  //                         arcs[a].end_point_id() ==
  //                             reinterpret_cast<std::uintptr_t>(&test_pt))
  //                             {
  //                       sprintf(faces_string[f], "%s%s -- ",
  //                       faces_string[f],
  //                               arcs_string[a]);
  //                       sprintf(faces_arrow_string[f],
  //                               "%s\\draw[very thick,-{angle "
  //                               "60[scale=0.75]},%s] %s;\n",
  //                               faces_arrow_string[f],
  //                               face_normal[2] > 0.0 ? "Red3" :
  //                               "Blue4", arcs_arrow_string[a]);
  //                       break;
  //                     } else if (arcs[a].end_point_id() ==
  //                                    reinterpret_cast<std::uintptr_t>(
  //                                        &next_pt) &&
  //                                arcs[a].start_point_id() ==
  //                                    reinterpret_cast<std::uintptr_t>(
  //                                        &test_pt)) {
  //                       sprintf(faces_string[f], "%s%s -- ",
  //                       faces_string[f],
  //                               arcs_string[nArcs + a]);
  //                       sprintf(faces_arrow_string[f],
  //                               "%s\\draw[very thick,-{angle "
  //                               "60[scale=0.75]},%s] %s;\n",
  //                               faces_arrow_string[f],
  //                               face_normal[2] > 0.0 ? "Red3" :
  //                               "Blue4",
  //                               arcs_arrow_string[nArcs + a]);
  //                       break;
  //                     }
  //                   }
  //                 }
  //               }
  //               if (current_half_edge == starting_half_edge) {
  //                 break;
  //               }
  //             }
  //           }
  //         }
  //         current_half_edge =
  //         current_half_edge->getNextHalfEdge();
  //       } while (current_half_edge != starting_half_edge);
  //       if (face_exists) {
  //         sprintf(faces_string[f], "%scycle", faces_string[f]);
  //       }
  //       center /= scale_center;
  //       std::cout << "\\clip " << faces_string[f] <<
  //       ";\n\\fill["
  //                 << ((face_normal[2] > 0.0) ? "DarkOrange4" :
  //                 "Blue4")
  //                 << "!30,opacity=0.4] " << faces_string[f] <<
  //                 ";\n\\draw["
  //                 << ((face_normal[2] > 0.0) ? "DarkOrange4" :
  //                 "Blue4")
  //                 << ",very thick] " << faces_string[f] <<
  //                 ";\n";
  //       // << faces_arrow_string[f]
  //       if (scale_center > 0.5) {
  //         std::cout << "\\draw["
  //                   << ((face_normal[2] > 0.0) ? "DarkOrange4"
  //                   : "Blue4")
  //                   << "] (" << center[0] << "," << center[1]
  //                   << ") node
  //                   {\\footnotesize$\\hat{\\mathcal{F}}^{\\perp}_{"
  //                   << f + 1 << "}$};\n";
  //       }
  //       std::cout << "\\end{scope}" << std::endl;
  //       std::cout << (face_normal[2] > 0.0 ?
  //       "\\end{pgfonlayer}{foreground}"
  //                                          :
  //                                          "\\end{pgfonlayer}{background}")
  //                 << std::endl;
  //     } else {
  //       unclipped_faces[f] = 0;
  //     }
  //   }
  //   std::cout << "%%%%%%%%%% END TIKZ CODE" << std::endl
  //             << std::endl
  //             << std::endl;
  // }
  // //////////////////////// FOR FIGURE

  return full_moments;
}
}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_TPP_
