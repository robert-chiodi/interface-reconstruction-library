// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_WEDGE_COMPUTATION_TPP_
#define IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_WEDGE_COMPUTATION_TPP_

#include <float.h>
#include <cassert>
#include <cmath>

#include <Eigen/Dense>

#include "irl/data_structures/stack_vector.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/pt.h"
#include "irl/helpers/mymath.h"
#include "irl/optimization/secant.h"
#include "irl/paraboloid_reconstruction/aligned_paraboloid.h"
#include "irl/paraboloid_reconstruction/parametrized_surface.h"

namespace IRL {

template <class ReturnType>
ReturnType calculateTriangleCorrection(const AlignedParaboloid& a_paraboloid,
                                       const Pt& a_pt_0, const Pt& a_pt_1,
                                       const Pt& a_pt_2);

template <>
inline Volume calculateTriangleCorrection(const AlignedParaboloid& a_paraboloid,
                                          const Pt& a_pt_0, const Pt& a_pt_1,
                                          const Pt& a_pt_2) {
  return (-a_paraboloid.a() * (a_pt_0[0] + a_pt_1[0]) *
              (a_pt_1[0] + a_pt_2[0]) +
          -a_paraboloid.b() * (a_pt_0[1] + a_pt_1[1]) *
              (a_pt_1[1] + a_pt_2[1]) -
          a_pt_0[2] - 2.0 * a_pt_1[2] - a_pt_2[2]) /
         12.0 *
         ((a_pt_1[1] - a_pt_2[1]) * a_pt_0[0] +
          (a_pt_2[1] - a_pt_0[1]) * a_pt_1[0] +
          (a_pt_0[1] - a_pt_1[1]) * a_pt_2[0]);
}

template <>
inline VolumeMoments calculateTriangleCorrection(
    const AlignedParaboloid& a_paraboloid, const Pt& a_pt_0, const Pt& a_pt_1,
    const Pt& a_pt_2) {
  auto moments = VolumeMoments::fromScalarConstant(0.0);
  std::cout << "Not yet implemented" << std::endl;
  std::exit(-1);
  return moments;
}

inline double signedDistance(const Pt& a_pt,
                             const AlignedParaboloid& a_paraboloid) {
  return a_paraboloid.a() * a_pt[0] * a_pt[0] +
         a_paraboloid.b() * a_pt[1] * a_pt[1] + a_pt[2];
}

inline Normal computeTangentVectorAtPoint(const AlignedParaboloid& a_paraboloid,
                                          const Plane& a_plane,
                                          const Pt& a_pt) {
  const auto gradF = getParaboloidSurfaceNormal(a_paraboloid, a_pt);
  Normal tangent_at_pt = crossProduct(a_plane.normal(), gradF);
  tangent_at_pt.normalize();
  return tangent_at_pt;
}

inline Normal computeAndCorrectTangentVectorAtPt(
    const AlignedParaboloid& a_paraboloid, const Plane& a_plane,
    const Pt& a_origin_pt, const Pt& a_end_pt, const Normal& a_end_tangent,
    const Pt& a_intersection_pt) {
  Normal tangent =
      computeTangentVectorAtPoint(a_paraboloid, a_plane, a_intersection_pt);
  const Normal correct_sign =
      crossProduct(a_plane.normal(), a_end_pt - a_origin_pt);
  if ((a_end_tangent * correct_sign > 0.0) !=
      (tangent * Normal(crossProduct(a_plane.normal(),
                                     a_intersection_pt - a_origin_pt)) >
       0.0)) {
    tangent = -tangent;
  }
  return tangent;
}

template <class HalfEdgeType>
bool orientInitialTangents(const HalfEdgeType* a_starting_half_edge,
                           const Plane& a_face_plane, Normal* a_tangent) {
  const auto& face_normal = a_face_plane.normal();
  Normal edge = a_starting_half_edge->getVertex()->getLocation() -
                a_starting_half_edge->getPreviousVertex()->getLocation();
  Normal edge_normal = crossProduct(face_normal, edge);
  const double dot = edge_normal * (*a_tangent);
  if (std::fabs(dot) > DBL_EPSILON) {
    if (dot < 0.0) {
      (*a_tangent) *= -1.0;
    }
    // First half edge non-zero length and not
    // colinear with face tangent.
    return true;
  }

  // The intersection now necessarily coincides with one of the vertices of the
  // edge
  auto previous_half_edge = a_starting_half_edge;
  double edge_length =
      squaredMagnitude(previous_half_edge->getVertex()->getLocation() -
                       previous_half_edge->getPreviousVertex()->getLocation());
  while (edge_length < MINIMUM_EDGE_LENGTH * MINIMUM_EDGE_LENGTH) {
    previous_half_edge = previous_half_edge->getPreviousHalfEdge();
    if (previous_half_edge == a_starting_half_edge) {
      return false;
    }
    edge_length = squaredMagnitude(
        previous_half_edge->getVertex()->getLocation() -
        previous_half_edge->getPreviousVertex()->getLocation());
  }

  auto next_half_edge = a_starting_half_edge->getNextHalfEdge();
  edge_length =
      squaredMagnitude(next_half_edge->getVertex()->getLocation() -
                       next_half_edge->getPreviousVertex()->getLocation());
  while (edge_length < MINIMUM_EDGE_LENGTH * MINIMUM_EDGE_LENGTH) {
    next_half_edge = next_half_edge->getNextHalfEdge();
    edge_length =
        squaredMagnitude(next_half_edge->getVertex()->getLocation() -
                         next_half_edge->getPreviousVertex()->getLocation());
  }

  edge = previous_half_edge->getVertex()->getLocation() -
         previous_half_edge->getPreviousVertex()->getLocation();
  edge_normal = crossProduct(face_normal, edge);
  double previous_dot_tangent = edge_normal * (*a_tangent);
  edge = next_half_edge->getVertex()->getLocation() -
         next_half_edge->getPreviousVertex()->getLocation();
  edge_normal = crossProduct(face_normal, edge);
  double next_dot_tangent = edge_normal * (*a_tangent);

  if (std::fabs(previous_dot_tangent) > DBL_EPSILON &&
      std::fabs(next_dot_tangent) > DBL_EPSILON) {
    if (previous_dot_tangent * next_dot_tangent > 0.0) {
      if (previous_dot_tangent < 0.0) {
        (*a_tangent) *= -1.0;
      }
      // Case where previous or next half edge was zero length, so had
      // to continue in both directions to get non-zero length edge
      // to orient
      return true;
    } else {
      // Edge normal vectors dotted with tangent have different sign,
      // meaning paraboloid intersects obliquely and no portion of
      // paraboloid actually lays inside face. Therefore, return
      // false and do not compute wedge correction.
      return false;
    }
  } else {
    edge = previous_half_edge->getPreviousVertex()->getLocation() -
           previous_half_edge->getVertex()->getLocation();
    previous_dot_tangent = edge * (*a_tangent);
    if (std::fabs(previous_dot_tangent) < DBL_EPSILON) {
      edge = next_half_edge->getVertex()->getLocation() -
             next_half_edge->getPreviousVertex()->getLocation();
    }
    if (edge * (*a_tangent) < 0.0) {
      (*a_tangent) *= -1.0;
    }
    // Case where Paraboloid intersects on a vertex and tangent is
    // colinear to one of the edges on that vertex.
    return true;
  }
}

inline std::array<double, 3> coeffsV3SeriesOne(const double a_weight) {
  auto coeffs = std::array<double, 3>({0.0, 0.0, 0.0});
  double x = 1.0;
  for (UnsignedIndex_t i = 0; i <= 40; ++i) {
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      coeffs[j] += v3Series[i][j] * x;
    }
    x *= a_weight - 1.0;
  }
  return coeffs;
}

inline std::array<double, 3> coeffsV3SeriesInfinity(const double a_weight) {
  const double wm2 = 1.0 / (a_weight * a_weight);
  const double wm4 = wm2 * wm2;
  const double wm6 = wm4 * wm2;
  const double wm8 = wm4 * wm4;
  const double ln2plnw = log(2.0) + log(a_weight);
  return std::array<double, 3>(
      {1.0 / 6.0 + wm2 / 4.0 + wm4 * (17.0 / 6.0 - ln2plnw * 7.0 / 4.0) +
           wm6 * (401.0 / 48.0 - ln2plnw * 55.0 / 8.0) +
           wm8 * (2225.0 / 128.0 - ln2plnw * 525.0 / 32.0),
       2.0 / 3.0 - wm2 + wm4 * (-23.0 / 3.0 + 5.0 * ln2plnw) +
           wm6 * (-247.0 / 12.0 + 35.0 * ln2plnw / 2.0) +
           wm8 * (-1307.0 / 32.0 + 315.0 * ln2plnw / 8.0),
       wm2 * (11.0 / 6.0 - ln2plnw) + wm4 * (77.0 / 12.0 - 5.0 * ln2plnw) +
           wm6 * (459.0 / 32.0 - 105.0 * ln2plnw / 8.0) +
           wm8 * (2509.0 / 96.0 - 105.0 * ln2plnw / 4.0)});
}

inline std::array<double, 3> coeffsV3Exact(const double a_weight) {
  const double w2 = a_weight * a_weight;
  const double w3 = w2 * a_weight;
  const double w4 = w2 * w2;
  const double w5 = w2 * w3;
  const double w6 = w3 * w3;
  const double L = 1.0 / ((a_weight - 1.0) * (a_weight + 1.0));
  const double L3 = L * L * L;
  const double S = (a_weight < 1.0) ? sqrt(1.0 - a_weight * a_weight)
                                    : sqrt(a_weight * a_weight - 1.0);
  const double T = (a_weight < 1.0) ? atan((1.0 - a_weight) / S) / S
                                    : atanh((a_weight - 1.0) / S) / S;
  return std::array<double, 3>(
      {(2.0 * w6 - 3.0 * w4 + 31.0 * w2 - (42.0 * w3 + 18.0 * a_weight) * T) *
           L3 / 12.0,
       (2.0 * w6 - 9.0 * w4 - 8.0 * w2 + 30.0 * w3 * T) * L3 / 3.0,
       (11.0 * w4 + 4.0 * w2 - (12.0 * w5 + 18.0 * w3) * T) * L3 / 6.0});
}

template <class ReturnType>
ReturnType computeV3Contribution(const AlignedParaboloid& a_paraboloid,
                                 const RationalBezierArc& a_arc);

template <>
inline Volume computeV3Contribution<Volume>(
    const AlignedParaboloid& a_paraboloid, const RationalBezierArc& a_arc) {
  const auto& pt_0 = a_arc.start_point();
  const auto& cp = a_arc.control_point();
  const auto& pt_1 = a_arc.end_point();
  const auto& weight = a_arc.weight();
  const double area_proj_triangle =
      0.5 * (pt_0[0] * (pt_1[1] - cp[1]) + pt_1[0] * (cp[1] - pt_0[1]) +
             cp[0] * (pt_0[1] - pt_1[1]));
  assert(weight >= 0.0);
  std::array<double, 3> coeffs;
  if (weight < 0.35)  // We use the exact expressions
    coeffs = coeffsV3Exact(weight);
  else if (weight < 1.7)  // We use the 40th order Taylor series (w -> 1)
    coeffs = coeffsV3SeriesOne(weight);
  else if (weight < 80.0)  // We use the exact expressions
    coeffs = coeffsV3Exact(weight);
  else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
    coeffs = coeffsV3SeriesInfinity(weight);
  else  // This is within DBL_EPSILON of the actual value
    coeffs = std::array<double, 3>({1.0 / 6.0, 2.0 / 3.0, 0.0});
  return area_proj_triangle *
         (coeffs[0] * signedDistance(0.5 * (pt_0 + pt_1), a_paraboloid) +
          coeffs[1] *
              signedDistance(0.25 * (pt_0 + pt_1) + 0.5 * cp, a_paraboloid) +
          coeffs[2] * signedDistance(cp, a_paraboloid));
}

template <>
inline VolumeMoments computeV3Contribution<VolumeMoments>(
    const AlignedParaboloid& a_paraboloid, const RationalBezierArc& a_arc) {
  auto moments = VolumeMoments::fromScalarConstant(0.0);
  std::cout << "Not yet implemented" << std::endl;
  std::exit(-1);
  return moments;
}

template <class ReturnType, class SurfaceOutputType>
ReturnType computeV3ContributionWithSplit(const AlignedParaboloid& a_paraboloid,
                                          const Plane& a_plane,
                                          const Pt& a_pt_ref, const Pt& a_pt_0,
                                          const Pt& a_pt_1,
                                          const Normal& a_tangent_0,
                                          const Normal& a_tangent_1,
                                          SurfaceOutputType* a_surface) {
  if (squaredMagnitude(a_pt_1 - a_pt_0) < 100.0 * DBL_EPSILON * DBL_EPSILON) {
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      a_surface->addArc(
          RationalBezierArc(a_pt_0, 0.5 * (a_pt_0 + a_pt_1), a_pt_1, 0.0));
    }
    return ReturnType::fromScalarConstant(0.0);
  } else {
    const Normal edge_vector = a_pt_1 - a_pt_0;
    bool split = ((a_tangent_0 * edge_vector) < 0.0 ||
                  (a_tangent_1 * edge_vector) > 0.0);
    if (split) {
      const Pt average_pt = 0.5 * (a_pt_0 + a_pt_1);
      auto average_tangent = Normal(0.5 * (a_tangent_0 + a_tangent_1));
      if (squaredMagnitude(average_tangent) <
          100.0 * DBL_EPSILON * DBL_EPSILON) {
        average_tangent = Normal(0.25 * a_tangent_0 + 0.75 * a_tangent_1);
      }
      average_tangent.normalize();
      Pt projected_pt = projectPtAlongHalfLineOntoParaboloid(
          a_paraboloid, average_tangent, average_pt);
      const Normal tangent_projected_pt = computeAndCorrectTangentVectorAtPt(
          a_paraboloid, a_plane, a_pt_0, a_pt_1, a_tangent_1, projected_pt);
      // We need to store this vertex so that its address remains unique over
      // time
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        Pt* new_point = new Pt(projected_pt);
        a_surface->addPt(new_point);
        return computeV3ContributionWithSplit<ReturnType>(
                   a_paraboloid, a_plane, a_pt_ref, a_pt_0, *new_point,
                   a_tangent_0, tangent_projected_pt, a_surface) +
               computeV3ContributionWithSplit<ReturnType>(
                   a_paraboloid, a_plane, a_pt_ref, *new_point, a_pt_1,
                   -tangent_projected_pt, a_tangent_1, a_surface);
      } else {
        return computeV3ContributionWithSplit<ReturnType>(
                   a_paraboloid, a_plane, a_pt_ref, a_pt_0, projected_pt,
                   a_tangent_0, tangent_projected_pt, a_surface) +
               computeV3ContributionWithSplit<ReturnType>(
                   a_paraboloid, a_plane, a_pt_ref, projected_pt, a_pt_1,
                   -tangent_projected_pt, a_tangent_1, a_surface);
      }
    } else {
      const auto arc = RationalBezierArc(a_pt_0, a_tangent_0, a_pt_1,
                                         a_tangent_1, a_plane, a_paraboloid);
      // Return surface parametrization
      if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        a_surface->addArc(arc);
      }
      return computeV3Contribution<ReturnType>(a_paraboloid, arc) +
             calculateTriangleCorrection<ReturnType>(a_paraboloid, a_pt_0,
                                                     a_pt_1, a_pt_ref);
    }
  }
}

template <class ReturnType, class SurfaceOutputType>
ReturnType bezierIntegrate(const AlignedParaboloid& a_paraboloid,
                           const Plane& a_plane, const Pt& a_pt_ref,
                           const Pt& a_pt_0, const Pt& a_pt_1,
                           const Normal& a_tangent_0, const Normal& a_tangent_1,
                           SurfaceOutputType* a_surface) {
  const Normal edge_vector = a_pt_1 - a_pt_0;
  const double dot_0 = a_tangent_0 * edge_vector;
  const double dot_1 = a_tangent_1 * edge_vector;
  const Pt average_pt = 0.5 * (a_pt_0 + a_pt_1);
  std::cout << dot_0 << " " << dot_1 << std::endl;
  std::cout << a_tangent_0 << " " << a_tangent_1 << std::endl;
  std::cout << a_pt_1 << " " << a_pt_0 << std::endl;
  std::cout << squaredMagnitude(a_pt_1 - a_pt_0) << std::endl;
  std::cout << a_plane.normal() << std::endl;
  if (squaredMagnitude(a_pt_1 - a_pt_0) < 100.0 * DBL_EPSILON * DBL_EPSILON) {
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      a_surface->addArc(RationalBezierArc(a_pt_1, average_pt, a_pt_0, 1.0));
    }
    return ReturnType::fromScalarConstant(0.0);
  } else if (std::fabs(a_tangent_0 * a_tangent_1 + 1.0) < 10.0 * DBL_EPSILON) {
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      a_surface->addArc(RationalBezierArc(a_pt_1, average_pt, a_pt_0, 1.0));
    }
    return ReturnType::fromScalarConstant(0.0);
  } else if (dot_0 <= 0.0 || dot_1 >= 0.0) {
    auto average_tangent = Normal(0.5 * (a_tangent_0 + a_tangent_1));
    Pt projected_pt = projectPtAlongHalfLineOntoParaboloid(
        a_paraboloid, average_tangent, average_pt);
    const Normal tangent_projected_pt = computeAndCorrectTangentVectorAtPt(
        a_paraboloid, a_plane, a_pt_0, a_pt_1, a_tangent_1, projected_pt);
    // We need to store this vertex so that its address remains unique over time
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      Pt* new_point = new Pt(projected_pt);
      a_surface->addPt(new_point);
      return bezierIntegrate<ReturnType>(a_paraboloid, a_plane, a_pt_ref,
                                         a_pt_0, *new_point, a_tangent_0,
                                         tangent_projected_pt, a_surface) +
             bezierIntegrate<ReturnType>(
                 a_paraboloid, a_plane, a_pt_ref, *new_point, a_pt_1,
                 -tangent_projected_pt, a_tangent_1, a_surface);
    } else {
      return bezierIntegrate<ReturnType>(a_paraboloid, a_plane, a_pt_ref,
                                         a_pt_0, projected_pt, a_tangent_0,
                                         tangent_projected_pt, a_surface) +
             bezierIntegrate<ReturnType>(
                 a_paraboloid, a_plane, a_pt_ref, projected_pt, a_pt_1,
                 -tangent_projected_pt, a_tangent_1, a_surface);
    }
  } else {
    const Normal n_cross_t0 = crossProduct(a_plane.normal(), a_tangent_0);
    assert(std::fabs(n_cross_t0 * a_tangent_1) > DBL_EPSILON);
    const double lambda_1 =
        -(n_cross_t0 * edge_vector) / (n_cross_t0 * a_tangent_1);
    const auto control_pt = Pt(a_pt_1 + lambda_1 * a_tangent_1);
    const auto mid_to_control = Normal(control_pt - average_pt);
    const auto projected_pt = projectPtAlongHalfLineOntoParaboloid(
        a_paraboloid, mid_to_control, average_pt);
    double weight = 1.0;
    if (squaredMagnitude(projected_pt - control_pt) <
        DBL_EPSILON * DBL_EPSILON) {
      weight = DBL_MAX;
    } else {
      double previous_best = -DBL_MAX;
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        const double denominator = projected_pt[d] - control_pt[d];
        if (std::fabs(denominator) > previous_best) {
          previous_best = std::fabs(denominator);
          weight = (a_pt_0[d] + a_pt_1[d] - 2.0 * projected_pt[d]) /
                   (2.0 * denominator);
        }
      }
    }
    // Return surface parametrization
    if constexpr (!std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
      a_surface->addArc(RationalBezierArc(a_pt_1, control_pt, a_pt_0, weight));
    }
    return computeV3Contribution<ReturnType>(
               a_paraboloid,
               RationalBezierArc(a_pt_0, control_pt, a_pt_1, weight)) +
           calculateTriangleCorrection<ReturnType>(a_paraboloid, a_pt_0, a_pt_1,
                                                   a_pt_ref);
  }
}

template <class ReturnType, class SurfaceOutputType>
inline ReturnType computeWedgeCorrection(const AlignedParaboloid& a_paraboloid,
                                         const Plane& a_plane, const Pt& a_pt_0,
                                         const Pt& a_pt_1,
                                         const Normal& a_tangent_0,
                                         const Normal& a_tangent_1,
                                         SurfaceOutputType* a_surface) {
  return bezierIntegrate<ReturnType>(a_paraboloid, a_plane, a_pt_0, a_pt_0,
                                     a_pt_1, a_tangent_0, a_tangent_1,
                                     a_surface);
}

}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_TPP_
