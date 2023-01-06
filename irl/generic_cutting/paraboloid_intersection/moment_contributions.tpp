// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_MOMENT_CONTRIBUTIONS_TPP_
#define IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_MOMENT_CONTRIBUTIONS_TPP_

#include <float.h>
#include <cassert>
#include <cmath>

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/scalar_with_gradient.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/helpers/mymath.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

namespace IRL {

/* This returns the algebraic signed distance to a paraboloid
 * returns < 0 number if a_pt is below the paraboloid
 * returns > 0 number if a_pt is above the paraboloid */
inline double signedDistance(const Pt& a_pt,
                             const AlignedParaboloid& a_paraboloid) {
  return a_paraboloid.a() * a_pt[0] * a_pt[0] +
         a_paraboloid.b() * a_pt[1] * a_pt[1] + a_pt[2];
}

/******************************************************************************/
/*********************** First moment contribution ****************************/
/******************************************************************************/
/* This compute the first contribution to the moments (arising from the
 * integration of the face plane primitives on the poligonized clipped faces) */
template <class ReturnType, class PtType>
ReturnType computeType1Contribution(const PtType& a_ref_pt,
                                    const PtType& a_pt_0, const PtType& a_pt_1);

/* This is the version of computeType1Contribution for returning the zeroth
 * moment only without its gradient */
template <>
inline enable_if_t<!has_embedded_gradient<Volume>::value, Volume>
computeType1Contribution(const Pt& a_ref_pt, const Pt& a_pt_0,
                         const Pt& a_pt_1) {
  return (a_ref_pt[2] + a_pt_0[2] + a_pt_1[2]) / 6.0 *
         ((a_pt_0[0] - a_ref_pt[0]) * (a_pt_1[1] - a_ref_pt[1]) -
          (a_pt_1[0] - a_ref_pt[0]) * (a_pt_0[1] - a_ref_pt[1]));
}

/* This is the version of computeType1Contribution for returning the zeroth
 * moment only or the zeroth and first moments with their gradient */
template <class ReturnType>
inline enable_if_t<has_embedded_gradient<ReturnType>::value, ReturnType>
computeType1Contribution(
    const PtWithGradient<typename ReturnType::gradient_type>& a_ref_pt,
    const PtWithGradient<typename ReturnType::gradient_type>& a_pt_0,
    const PtWithGradient<typename ReturnType::gradient_type>& a_pt_1) {
  using gradient_type = typename ReturnType::gradient_type;
  using scalar_type = ScalarWithGradient<gradient_type>;
  auto moments_with_gradient = ReturnType::fromScalarConstant(0.0);
  std::array<scalar_type, 3> pt_0, pt_1, ref_pt;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_0[d] = scalar_type(a_pt_0.getPt()[d], a_pt_0.getData()[d]);
    pt_1[d] = scalar_type(a_pt_1.getPt()[d], a_pt_1.getData()[d]);
    ref_pt[d] = scalar_type(a_ref_pt.getPt()[d], a_ref_pt.getData()[d]);
  }
  const auto triangle_area = ((pt_0[0] - ref_pt[0]) * (pt_1[1] - ref_pt[1]) -
                              (pt_1[0] - ref_pt[0]) * (pt_0[1] - ref_pt[1])) /
                             2.0;
  const auto volume = triangle_area * (ref_pt[2] + pt_0[2] + pt_1[2]) / 3.0;
  moments_with_gradient.volume() = volume.value();
  moments_with_gradient.volume_gradient() = volume.gradient();
  if constexpr (std::is_same<typename ReturnType::moment_type,
                             VolumeMoments>::value) {
    Pt& centroid = moments_with_gradient.centroid().getPt();
    auto& centroid_gradient = moments_with_gradient.centroid().getData();
    const auto m1x = triangle_area *
                     (pt_0[0] * (2.0 * pt_0[2] + pt_1[2] + ref_pt[2]) +
                      pt_1[0] * (pt_0[2] + 2.0 * pt_1[2] + ref_pt[2]) +
                      ref_pt[0] * (pt_0[2] + pt_1[2] + 2.0 * ref_pt[2])) /
                     12.0;
    const auto m1y = triangle_area *
                     (pt_0[1] * (2.0 * pt_0[2] + pt_1[2] + ref_pt[2]) +
                      pt_1[1] * (pt_0[2] + 2.0 * pt_1[2] + ref_pt[2]) +
                      ref_pt[1] * (pt_0[2] + pt_1[2] + 2.0 * ref_pt[2])) /
                     12.0;
    const auto m1z =
        triangle_area *
        (pt_0[2] * pt_0[2] + pt_1[2] * pt_1[2] + ref_pt[2] * ref_pt[2] +
         pt_1[2] * ref_pt[2] + pt_0[2] * pt_1[2] + pt_0[2] * ref_pt[2]) /
        12.0;
    centroid[0] = m1x.value();
    centroid[1] = m1y.value();
    centroid[2] = m1z.value();
    centroid_gradient[0] = m1x.gradient();
    centroid_gradient[1] = m1y.gradient();
    centroid_gradient[2] = m1z.gradient();
  }
  return moments_with_gradient;
}

/* This is the version of computeType1Contribution for returning the zeroth and
 * first moments without their gradient */
template <>
inline enable_if_t<!has_embedded_gradient<VolumeMoments>::value, VolumeMoments>
computeType1Contribution(const Pt& a_ref_pt, const Pt& a_pt_0,
                         const Pt& a_pt_1) {
  auto moments = VolumeMoments::fromScalarConstant(0.0);
  const double triangle_area =
      ((a_pt_0[0] - a_ref_pt[0]) * (a_pt_1[1] - a_ref_pt[1]) -
       (a_pt_1[0] - a_ref_pt[0]) * (a_pt_0[1] - a_ref_pt[1])) /
      2.0;
  moments.volume() =
      triangle_area * (a_ref_pt[2] + a_pt_0[2] + a_pt_1[2]) / 3.0;
  moments.centroid()[0] =
      triangle_area *
      (a_pt_0[0] * (2.0 * a_pt_0[2] + a_pt_1[2] + a_ref_pt[2]) +
       a_pt_1[0] * (a_pt_0[2] + 2.0 * a_pt_1[2] + a_ref_pt[2]) +
       a_ref_pt[0] * (a_pt_0[2] + a_pt_1[2] + 2.0 * a_ref_pt[2])) /
      12.0;
  moments.centroid()[1] =
      triangle_area *
      (a_pt_0[1] * (2.0 * a_pt_0[2] + a_pt_1[2] + a_ref_pt[2]) +
       a_pt_1[1] * (a_pt_0[2] + 2.0 * a_pt_1[2] + a_ref_pt[2]) +
       a_ref_pt[1] * (a_pt_0[2] + a_pt_1[2] + 2.0 * a_ref_pt[2])) /
      12.0;
  moments.centroid()[2] = triangle_area *
                          (a_pt_0[2] * a_pt_0[2] + a_pt_1[2] * a_pt_1[2] +
                           a_ref_pt[2] * a_ref_pt[2] + a_pt_1[2] * a_ref_pt[2] +
                           a_pt_0[2] * a_pt_1[2] + a_pt_0[2] * a_ref_pt[2]) /
                          12.0;
  return moments;
}

template <class ReturnType, class PtType>
ReturnType computeType2Contribution(
    const AlignedParaboloid& a_aligned_paraboloid, const PtType& a_pt_0,
    const PtType& a_pt_1);

template <>
inline enable_if_t<!has_embedded_gradient<Volume>::value, Volume>
computeType2Contribution(const AlignedParaboloid& a_aligned_paraboloid,
                         const Pt& a_pt_0, const Pt& a_pt_1) {
  return (a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) / 12.0 *
         (-a_pt_0[2] - a_pt_1[2] +
          a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] +
          a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1]);
}

template <class ReturnType>
inline enable_if_t<has_embedded_gradient<ReturnType>::value, ReturnType>
computeType2Contribution(
    const AlignedParaboloid& a_aligned_paraboloid,
    const PtWithGradient<typename ReturnType::gradient_type>& a_pt_0,
    const PtWithGradient<typename ReturnType::gradient_type>& a_pt_1) {
  using gradient_type = typename ReturnType::gradient_type;
  using scalar_type = ScalarWithGradient<gradient_type>;
  auto moments_with_gradient = ReturnType::fromScalarConstant(0.0);
  std::array<scalar_type, 3> pt_0, pt_1;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_0[d] = scalar_type(a_pt_0.getPt()[d], a_pt_0.getData()[d]);
    pt_1[d] = scalar_type(a_pt_1.getPt()[d], a_pt_1.getData()[d]);
  }
  auto A_grad = gradient_type(0.0), B_grad = gradient_type(0.0);
  A_grad.setGradA(1.0);
  B_grad.setGradB(1.0);
  const auto A = scalar_type(a_aligned_paraboloid.a(), A_grad),
             B = scalar_type(a_aligned_paraboloid.b(), B_grad);
  const auto volume =
      (pt_0[0] * pt_1[1] - pt_1[0] * pt_0[1]) / 12.0 *
      (-pt_0[2] - pt_1[2] + A * pt_0[0] * pt_1[0] + B * pt_0[1] * pt_1[1]);
  moments_with_gradient.volume() = volume.value();
  moments_with_gradient.volume_gradient() = volume.gradient();
  if constexpr (std::is_same<typename ReturnType::moment_type,
                             VolumeMoments>::value) {
    Pt& centroid = moments_with_gradient.centroid().getPt();
    auto& centroid_gradient = moments_with_gradient.centroid().getData();
    const auto m1x = (pt_1[0] * pt_0[1] - pt_0[0] * pt_1[1]) *
                     (2.0 * B * (pt_0[1] - pt_1[1]) *
                          (pt_1[0] * pt_0[1] - pt_0[0] * pt_1[1]) +
                      3.0 * (pt_0[0] + pt_1[0]) * (pt_0[2] + pt_1[2])) /
                     60.0;
    const auto m1y = (pt_1[0] * pt_0[1] - pt_0[0] * pt_1[1]) *
                     (2.0 * A * (pt_0[0] - pt_1[0]) *
                          (pt_1[1] * pt_0[0] - pt_0[1] * pt_1[0]) +
                      3.0 * (pt_0[1] + pt_1[1]) * (pt_0[2] + pt_1[2])) /
                     60.0;
    const auto m1z =
        ((pt_0[0] * pt_1[1] - pt_1[0] * pt_0[1]) *
         (2.0 * A * B *
              ((pt_1[0] * pt_0[1] - pt_0[0] * pt_1[1]) *
               (pt_1[0] * pt_0[1] - pt_0[0] * pt_1[1])) +
          3.0 * A * pt_0[0] * pt_1[0] * (pt_0[2] + pt_1[2]) +
          3.0 * B * pt_0[1] * pt_1[1] * (pt_0[2] + pt_1[2]) -
          3.0 * (pt_0[2] * pt_0[2] + pt_0[2] * pt_1[2] + pt_1[2] * pt_1[2]))) /
        180.0;
    centroid[0] = m1x.value();
    centroid[1] = m1y.value();
    centroid[2] = m1z.value();
    centroid_gradient[0] = m1x.gradient();
    centroid_gradient[1] = m1y.gradient();
    centroid_gradient[2] = m1z.gradient();
  }
  return moments_with_gradient;
}

template <>
inline enable_if_t<!has_embedded_gradient<VolumeMoments>::value, VolumeMoments>
computeType2Contribution(const AlignedParaboloid& a_aligned_paraboloid,
                         const Pt& a_pt_0, const Pt& a_pt_1) {
  auto moments = VolumeMoments::fromScalarConstant(0.0);
  moments.volume() = (a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) / 12.0 *
                     (-a_pt_0[2] - a_pt_1[2] +
                      a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] +
                      a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1]);
  moments.centroid()[0] =
      (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) *
      (2.0 * a_aligned_paraboloid.b() * (a_pt_0[1] - a_pt_1[1]) *
           (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) +
       3.0 * (a_pt_0[0] + a_pt_1[0]) * (a_pt_0[2] + a_pt_1[2])) /
      60.0;
  moments.centroid()[1] =
      (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) *
      (2.0 * a_aligned_paraboloid.a() * (a_pt_0[0] - a_pt_1[0]) *
           (a_pt_1[1] * a_pt_0[0] - a_pt_0[1] * a_pt_1[0]) +
       3.0 * (a_pt_0[1] + a_pt_1[1]) * (a_pt_0[2] + a_pt_1[2])) /
      60.0;
  moments.centroid()[2] =
      ((a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) *
       (2.0 * a_aligned_paraboloid.a() * a_aligned_paraboloid.b() *
            ((a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) *
             (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1])) +
        3.0 * a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] *
            (a_pt_0[2] + a_pt_1[2]) +
        3.0 * a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1] *
            (a_pt_0[2] + a_pt_1[2]) -
        3.0 * (a_pt_0[2] * a_pt_0[2] + a_pt_0[2] * a_pt_1[2] +
               a_pt_1[2] * a_pt_1[2]))) /
      180.0;
  return moments;
}

template <class ScalarType>
inline std::array<ScalarType, 3> coeffsV3SeriesOne(const ScalarType& a_weight) {
  std::array<ScalarType, 3> coeffs;
  coeffs.fill(ScalarType(0.0));
  ScalarType x(1.0);
  for (UnsignedIndex_t i = 0; i <= 40; ++i) {
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      coeffs[j] += v3Series[i][j] * x;
    }
    x *= a_weight - ScalarType(1.0);
  }
  return coeffs;
}

template <class ScalarType>
inline std::array<ScalarType, 12> coeffsV3andC3SeriesOne(
    const ScalarType& a_weight) {
  std::array<ScalarType, 12> coeffs;
  coeffs.fill(ScalarType(0.0));
  ScalarType x(1.0);
  for (UnsignedIndex_t i = 0; i <= 40; ++i) {
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      coeffs[j] += v3Series[i][j] * x;
    }
    for (UnsignedIndex_t j = 0; j < 4; ++j) {
      coeffs[3 + j] += cx3Series[i][j] * x;
    }
    for (UnsignedIndex_t j = 0; j < 5; ++j) {
      coeffs[7 + j] += cz3Series[i][j] * x;
    }
    x *= a_weight - ScalarType(1.0);
  }
  return coeffs;
}

template <class ScalarType>
inline std::array<ScalarType, 3> coeffsV3SeriesInfinity(
    const ScalarType& a_weight) {
  const auto wm2 = 1.0 / (a_weight * a_weight);
  const auto wm4 = wm2 * wm2;
  const auto wm6 = wm4 * wm2;
  const auto wm8 = wm4 * wm4;
  const auto ln2plnw =
      ScalarType(std::log(2.0)) + LogMoments<ScalarType>(a_weight);
  return std::array<ScalarType, 3>(
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

template <class ScalarType>
inline std::array<ScalarType, 12> coeffsV3andC3SeriesInfinity(
    const double a_weight) {
  std::cout << "Not implemented yet" << std::endl;
  std::exit(-1);
  return std::array<ScalarType, 12>();
}

template <class ScalarType>
inline std::array<ScalarType, 3> coeffsV3Exact(const ScalarType& a_weight) {
  const auto w2 = a_weight * a_weight;
  const auto w3 = w2 * a_weight;
  const auto w4 = w2 * w2;
  const auto w5 = w2 * w3;
  const auto w6 = w3 * w3;
  const auto L = ScalarType(1.0) /
                 ((a_weight - ScalarType(1.0)) * (a_weight + ScalarType(1.0)));
  const auto L3 = L * L * L;
  const auto S =
      (a_weight < ScalarType(1.0))
          ? SqrtMoments<ScalarType>(ScalarType(1.0) - a_weight * a_weight)
          : SqrtMoments<ScalarType>(a_weight * a_weight - ScalarType(1.0));
  const auto T =
      (a_weight < ScalarType(1.0))
          ? ArctanMoments<ScalarType>((ScalarType(1.0) - a_weight) / S) / S
          : ArctanhMoments<ScalarType>((a_weight - ScalarType(1.0)) / S) / S;
  return std::array<ScalarType, 3>(
      {(2.0 * w6 - 3.0 * w4 + 31.0 * w2 - (42.0 * w3 + 18.0 * a_weight) * T) *
           L3 / 12.0,
       (2.0 * w6 - 9.0 * w4 - 8.0 * w2 + 30.0 * w3 * T) * L3 / 3.0,
       (11.0 * w4 + 4.0 * w2 - (12.0 * w5 + 18.0 * w3) * T) * L3 / 6.0});
}

template <class ScalarType>
inline std::array<ScalarType, 12> coeffsV3andC3Exact(
    const ScalarType& a_weight) {
  const auto w2 = a_weight * a_weight;
  const auto w3 = w2 * a_weight;
  const auto w4 = w2 * w2;
  const auto w5 = w2 * w3;
  const auto w6 = w3 * w3;
  const auto w7 = w4 * w3;
  const auto w8 = w4 * w4;
  const auto w9 = w4 * w5;
  const auto w10 = w5 * w5;
  const auto L = ScalarType(1.0) /
                 ((a_weight - ScalarType(1.0)) * (a_weight + ScalarType(1.0)));
  const auto L3 = L * L * L;
  const auto L4 = L3 * L;
  const auto L5 = L4 * L;
  const auto S =
      (a_weight < ScalarType(1.0))
          ? SqrtMoments<ScalarType>(ScalarType(1.0) - a_weight * a_weight)
          : SqrtMoments<ScalarType>(a_weight * a_weight - ScalarType(1.0));
  const auto T =
      (a_weight < ScalarType(1.0))
          ? ArctanMoments<ScalarType>((ScalarType(1.0) - a_weight) / S) / S
          : ArctanhMoments<ScalarType>((a_weight - ScalarType(1.0)) / S) / S;
  return std::array<ScalarType, 12>(
      {L3 *
           (2.0 * w6 - 3.0 * w4 + 31.0 * w2 -
            (42.0 * w3 + 18.0 * a_weight) * T) /
           12.0,
       L3 * (2.0 * w6 - 9.0 * w4 - 8.0 * w2 + 30.0 * w3 * T) / 3.0,
       L3 * (11.0 * w4 + 4.0 * w2 - (12.0 * w5 + 18.0 * w3) * T) / 6.0,
       L4 * ((-T * a_weight) / 32. + (93. * (w2)) / 2240. -
             (163. * (w4)) / 3360. + (5. * (w6)) / 168. - (w8) / 140.),
       L4 * ((w2) / 70. + (-T * (w3)) / 16. + (29. * (w4)) / 1120. -
             (19. * (w6)) / 1680. + (w8) / 420.),
       -L4 * ((w2) / 210. - (w4) / 21. - (-T * (w5)) / 8. -
              (13. * (w6)) / 560. + (w8) / 280.),
       L4 * ((w2) / 35. - (16. * (w4)) / 105. + (58. * (w6)) / 105. - T * (w7) +
             (w8) / 14.),
       L5 * ((-T * a_weight) / 128. + (193. * (w2)) / 16128. -
             (149. * (w4)) / 8064. + (19. * (w6)) / 1120. -
             (41. * (w8)) / 5040. + (w10) / 630.),
       L5 * ((4. * (w2)) / 945. + (-T * (w3)) / 48. + (65. * (w4)) / 6048. -
             (w6) / 144. + (11. * (w8)) / 3780. - (w10) / 1890.),
       -L5 * ((w2) / 1890. - (13. * (w4)) / 1890. - (-T * (w5)) / 48. -
              (11. * (w6)) / 2016. + (5. * (w8)) / 3024. - (w10) / 3780.),
       L5 * ((w2) / 315. - (w4) / 45. + (4. * (w6)) / 35. + (-T * (w7)) / 4. +
             (17. * (w8)) / 504. - (w10) / 252.),
       -L5 * ((w2) / 63. - (29. * (w4)) / 315. + (26. * (w6)) / 105. -
              (194. * (w8)) / 315. + T * (w9) - (w10) / 18.)});
}

template <class ReturnType, class RationalBezierArcType>
ReturnType computeType3Contribution(const AlignedParaboloid& a_paraboloid,
                                    const RationalBezierArcType& a_arc);

template <>
inline enable_if_t<!has_embedded_gradient<Volume>::value, Volume>
computeType3Contribution(const AlignedParaboloid& a_paraboloid,
                         const RationalBezierArc& a_arc) {
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
    coeffs = coeffsV3Exact<double>(weight);
  else if (weight < 1.7)  // We use the 40th order Taylor series (w -> 1)
    coeffs = coeffsV3SeriesOne<double>(weight);
  else if (weight < 80.0)  // We use the exact expressions
    coeffs = coeffsV3Exact<double>(weight);
  else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
    coeffs = coeffsV3SeriesInfinity<double>(weight);
  else  // This is within DBL_EPSILON of the actual value
    coeffs = std::array<double, 3>({1.0 / 6.0, 2.0 / 3.0, 0.0});
  return area_proj_triangle *
         (coeffs[0] * signedDistance(0.5 * (pt_0 + pt_1), a_paraboloid) +
          coeffs[1] *
              signedDistance(0.25 * (pt_0 + pt_1) + 0.5 * cp, a_paraboloid) +
          coeffs[2] * signedDistance(cp, a_paraboloid));
}

template <class ReturnType>
inline enable_if_t<has_embedded_gradient<ReturnType>::value, ReturnType>
computeType3Contribution(
    const AlignedParaboloid& a_paraboloid,
    const RationalBezierArcWithGradient<
        PtWithGradient<typename ReturnType::gradient_type>>& a_arc) {
  using gradient_type = typename ReturnType::gradient_type;
  using scalar_type = ScalarWithGradient<gradient_type>;
  auto moments_with_gradient = ReturnType::fromScalarConstant(0.0);
  std::array<scalar_type, 3> pt_0, cp, pt_1;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_0[d] = scalar_type(a_arc.start_point().getPt()[d],
                          a_arc.start_point().getData()[d]);
    cp[d] = scalar_type(a_arc.control_point().getPt()[d],
                        a_arc.control_point().getData()[d]);
    pt_1[d] = scalar_type(a_arc.end_point().getPt()[d],
                          a_arc.end_point().getData()[d]);
  }
  auto A_grad = gradient_type(0.0), B_grad = gradient_type(0.0);
  A_grad.setGradA(1.0);
  B_grad.setGradB(1.0);
  const auto A = scalar_type(a_paraboloid.a(), A_grad),
             B = scalar_type(a_paraboloid.b(), B_grad);
  const auto X0 = pt_0[0], X1 = cp[0], X2 = pt_1[0];
  const auto Y0 = pt_0[1], Y1 = cp[1], Y2 = pt_1[1];
  const auto Z0 = pt_0[2], Z1p = cp[2], Z2 = pt_1[2];
  const auto area_proj_triangle =
      0.5 * (X0 * (Y2 - Y1) + X2 * (Y1 - Y0) + X1 * (Y0 - Y2));
  const auto m0_basis = std::array<scalar_type, 3>(
      {area_proj_triangle *
           (0.25 * A * (X0 + X2) * (X0 + X2) +
            0.25 * B * (Y0 + Y2) * (Y0 + Y2) + 0.5 * (Z0 + Z2)),
       area_proj_triangle *
           (A * (0.25 * (X0 + X2) + 0.5 * X1) * (0.25 * (X0 + X2) + 0.5 * X1) +
            B * (0.25 * (Y0 + Y2) + 0.5 * Y1) * (0.25 * (Y0 + Y2) + 0.5 * Y1) +
            (0.25 * (Z0 + Z2) + 0.5 * Z1p)),
       area_proj_triangle * (A * X1 * X1 + B * Y1 * Y1 + Z1p)});
  const auto weight = scalar_type(a_arc.weight(), a_arc.weight_gradient());
  assert(weight.value() >= 0.0);
  if constexpr (std::is_same<typename ReturnType::moment_type, Volume>::value) {
    std::array<scalar_type, 3> coeffs;
    if (weight < scalar_type(0.35))  // We use the exact expressions
      coeffs = coeffsV3Exact<scalar_type>(weight);
    else if (weight <
             scalar_type(1.7))  // We use the 40th order Taylor series (w -> 1)
      coeffs = coeffsV3SeriesOne<scalar_type>(weight);
    else if (weight <
             scalar_type(1.0e9))  // We use the series expansion (w -> infty)
      coeffs = coeffsV3Exact<scalar_type>(weight);
    else  // This is within DBL_EPSILON of the actual value
      coeffs = std::array<scalar_type, 3>(
          {scalar_type(1.0 / 6.0), scalar_type(2.0 / 3.0), scalar_type(0.0)});
    for (size_t i = 0; i < 3; ++i) {
      moments_with_gradient.volume() += coeffs[i].value() * m0_basis[i].value();
      moments_with_gradient.volume_gradient() +=
          coeffs[i].value() * m0_basis[i].gradient() +
          coeffs[i].gradient() * m0_basis[i].value();
    }
  } else {
    const auto Z1P = -A * X1 * X1 - B * Y1 * Y1;
    const auto AA = A * A, BB = B * B, AB = A * B;
    const auto X00 = X0 * X0, X11 = X1 * X1, X22 = X2 * X2;
    const auto X000 = X00 * X0, X111 = X11 * X1, X222 = X22 * X2;
    const auto X0000 = X00 * X00, X1111 = X11 * X11, X2222 = X22 * X22;
    const auto Y00 = Y0 * Y0, Y11 = Y1 * Y1, Y22 = Y2 * Y2;
    const auto Y000 = Y00 * Y0, Y111 = Y11 * Y1, Y222 = Y22 * Y2;
    const auto Y0000 = Y00 * Y00, Y1111 = Y11 * Y11, Y2222 = Y22 * Y22;
    const auto Z00 = Z0 * Z0, Z22 = Z2 * Z2;
    const auto Z1p1p = Z1p * Z1p, Z1P1P = Z1P * Z1P;
    const auto X02 = X0 * X2, X12 = X1 * X2, X01 = X0 * X1;
    const auto Y02 = Y0 * Y2, Y12 = Y1 * Y2, Y01 = Y0 * Y1;
    const auto Z02 = Z0 * Z2, Z01p = Z0 * Z1p, Z01P = Z0 * Z1P;
    const auto Z1p1P = Z1p * Z1P, Z1p2 = Z1p * Z2, Z1P2 = Z1P * Z2;
    const auto X0Z0 = X0 * Z0, X0Z1p = X0 * Z1p, X0Z1P = X0 * Z1P,
               X0Z2 = X0 * Z2;
    const auto X1Z0 = X1 * Z0, X1Z1p = X1 * Z1p, X1Z1P = X1 * Z1P,
               X1Z2 = X1 * Z2;
    const auto X2Z0 = X2 * Z0, X2Z1p = X2 * Z1p, X2Z1P = X2 * Z1P,
               X2Z2 = X2 * Z2;
    const auto Y0Z0 = Y0 * Z0, Y0Z1p = Y0 * Z1p, Y0Z1P = Y0 * Z1P,
               Y0Z2 = Y0 * Z2;
    const auto Y1Z0 = Y1 * Z0, Y1Z1p = Y1 * Z1p, Y1Z1P = Y1 * Z1P,
               Y1Z2 = Y1 * Z2;
    const auto Y2Z0 = Y2 * Z0, Y2Z1p = Y2 * Z1p, Y2Z1P = Y2 * Z1P,
               Y2Z2 = Y2 * Z2;
    // Compute coefficients (functions of the weight)
    std::array<scalar_type, 12> coeffs;
    if (weight < scalar_type(0.35))  // We use the exact expressions
    {
      coeffs = coeffsV3andC3Exact<scalar_type>(weight);
    } else if (weight < scalar_type(1.7))  // We use Taylor series (w -> 1)
    {
      coeffs = coeffsV3andC3SeriesOne<scalar_type>(weight);
    } else if (weight < scalar_type(1.0e9))  // We use the exact expressions
    {
      coeffs = coeffsV3andC3Exact<scalar_type>(weight);
    }
    // else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
    //   coeffs = coeffsV3andC3SeriesInfinity(weight);
    else  // This is within DBL_EPSILON of the actual value
    {
      coeffs = std::array<scalar_type, 12>(
          {scalar_type(1.0 / 6.0), scalar_type(2.0 / 3.0),
           scalar_type(-1.0 / 140.0), scalar_type(1.0 / 420.0),
           scalar_type(-1.0 / 280.0), scalar_type(1.0 / 14.0), scalar_type(0.0),
           scalar_type(1.0 / 630.0), scalar_type(-1.0 / 1890.0),
           scalar_type(1.0 / 3780.0), scalar_type(-1.0 / 252.0),
           scalar_type(1.0 / 18.0)});
    }
    auto m1x_basis = std::array<scalar_type, 4>(
        {-6. * (X0Z0 - X0Z2 - X2Z0 + X2Z2 - 2. * B * X2 * Y00 +
                2. * B * X0 * Y02 + 2. * B * X2 * Y02 - 2. * B * X0 * Y22),
         2. * (5. * X0Z0 + 10. * X0Z1p + 6. * X0Z1P + 7. * X0Z2 +
               30. * A * X02 * X1 - 11. * X1Z0 - 4. * X1Z1p - 11. * X1Z2 +
               7. * X2Z0 + 10. * X2Z1p + 6. * X2Z1P + 5. * X2Z2 -
               14. * B * X1 * Y00 + 4. * B * X2 * Y00 + 14. * B * X0 * Y01 -
               4. * B * X1 * Y01 + 10. * B * X2 * Y01 - 4. * B * X0 * Y02 +
               10. * B * X1 * Y02 - 4. * B * X2 * Y02 + 4. * B * X0 * Y11 +
               4. * B * X2 * Y11 + 10. * B * X0 * Y12 - 4. * B * X1 * Y12 +
               14. * B * X2 * Y12 + 4. * B * X0 * Y22 - 14. * B * X1 * Y22),
         2. * (-5. * X0Z1p + 18. * X0Z1P + X0Z2 + 6. * A * X02 * X1 -
               5. * X1Z0 - 6. * X1Z1p - 6. * X1Z1P - 5. * X1Z2 + X2Z0 -
               5. * X2Z1p + 18. * X2Z1P - 12. * B * X1 * Y01 +
               2. * B * X2 * Y01 + 2. * B * X1 * Y02 + 12. * B * X0 * Y11 +
               12. * B * X2 * Y11 + 2. * B * X0 * Y12 - 12. * B * X1 * Y12),
         2. * (X1Z1p - X1Z1P)});
    auto m1y_basis = std::array<scalar_type, 4>(
        {6. * (-Y0Z0 + Y0Z2 + 2. * A * (X22 * Y0 + X00 * Y2 - X02 * (Y0 + Y2)) +
               Y2Z0 - Y2Z2),
         2. * (5. * Y0Z0 + 10. * Y0Z1p + 6. * Y0Z1P + 7. * Y0Z2 +
               30. * B * Y02 * Y1 - 11. * Y1Z0 - 4. * Y1Z1p - 11. * Y1Z2 +
               2. * A *
                   (-2. * X02 * Y0 + 2. * X11 * Y0 + 5. * X12 * Y0 +
                    2. * X22 * Y0 - 7. * X00 * Y1 + 5. * X02 * Y1 -
                    2. * X12 * Y1 - 7. * X22 * Y1 + 2. * X00 * Y2 -
                    2. * X02 * Y2 + 2. * X11 * Y2 + 7. * X12 * Y2 +
                    X01 * (7. * Y0 - 2. * Y1 + 5. * Y2)) +
               7. * Y2Z0 + 10. * Y2Z1p + 6. * Y2Z1P + 5. * Y2Z2),
         -2. * (5. * Y0Z1p - 18. * Y0Z1P - Y0Z2 - 6. * B * Y02 * Y1 +
                5. * Y1Z0 + 6. * Y1Z1p + 6. * Y1Z1P + 5. * Y1Z2 -
                2. * A *
                    (X12 * Y0 - 6. * X01 * Y1 + X02 * Y1 - 6. * X12 * Y1 +
                     X01 * Y2 + 6. * X11 * (Y0 + Y2)) -
                Y2Z0 + 5. * Y2Z1p - 18. * Y2Z1P),
         2. * (Y1Z1p - Y1Z1P)});
    auto m1z_basis = std::array<scalar_type, 5>(
        {-(AA * (21. * X0000 + 28. * X000 * X2 + 30. * X00 * X22 +
                 28. * X0 * X222 + 21. * X2222)) -
             21. * BB * Y0000 - 28. * BB * Y000 * Y2 - 30. * BB * Y00 * Y22 -
             2. * AB *
                 (X00 * (21. * Y00 + 14. * Y02 + 5. * Y22) +
                  2. * X02 * (7. * Y00 + 10. * Y02 + 7. * Y22) +
                  X22 * (5. * Y00 + 14. * Y02 + 21. * Y22)) -
             28. * BB * Y0 * Y222 - 21. * BB * Y2222 + 40. * Z00 + 48. * Z02 +
             40. * Z22,
         3. * AA *
                 (21. * X000 * X1 - 7. * X00 * X11 - 10. * X02 * X11 +
                  35. * X00 * X12 - 7. * X000 * X2 - 10. * X00 * X22 +
                  35. * X01 * X22 - 7. * X11 * X22 - 7. * X0 * X222 +
                  21. * X1 * X222) -
             AB * (7. * X11 * Y00 - 35. * X12 * Y00 + 10. * X22 * Y00 -
                   63. * X00 * Y01 + 20. * X12 * Y01 - 35. * X22 * Y01 +
                   21. * X00 * Y02 + 10. * X11 * Y02 - 70. * X12 * Y02 +
                   21. * X22 * Y02 + 7. * X00 * Y11 + 7. * X22 * Y11 -
                   35. * X00 * Y12 + 28. * X12 * Y12 - 63. * X22 * Y12 +
                   X01 * (-63. * Y00 + 28. * Y01 - 70. * Y02 + 20. * Y12 -
                          35. * Y22) +
                   10. * X00 * Y22 + 7. * X11 * Y22 - 63. * X12 * Y22 +
                   X02 * (21. * Y00 - 70. * Y01 + 40. * Y02 + 10. * Y11 -
                          70. * Y12 + 21. * Y22)) -
             3. * (BB * (7. * Y00 * Y11 + 10. * Y02 * Y11 - 35. * Y00 * Y12 +
                         7. * Y000 * (-3. * Y1 + Y2) + 10. * Y00 * Y22 -
                         35. * Y01 * Y22 + 7. * Y11 * Y22 + 7. * Y0 * Y222 -
                         21. * Y1 * Y222) +
                   2. * (5. * Z00 + 10. * Z01p + 4. * Z02 - 2. * Z1p1p +
                         10. * Z1p2 + 5. * Z22)),
         -6. * AA *
                 (46. * X02 * X11 - 14. * X0 * X111 + X1111 - 14. * X111 * X2 -
                  14. * X01 * X22 + 28. * X11 * X22 +
                  X00 * (28. * X11 - 14. * X12 + X22)) -
             2. * AB *
                 (X22 * Y00 + 112. * X01 * Y01 - 28. * X02 * Y01 -
                  14. * X22 * Y01 - 28. * X01 * Y02 + 4. * X02 * Y02 +
                  28. * X00 * Y11 - 42. * X01 * Y11 + 46. * X02 * Y11 +
                  28. * X22 * Y11 -
                  2. * X12 *
                      (7. * Y00 - 46. * Y01 + 14. * Y02 + 21. * Y11 -
                       56. * Y12) -
                  14. * X00 * Y12 + 92. * X01 * Y12 - 28. * X02 * Y12 +
                  X00 * Y22 - 14. * X01 * Y22 +
                  X11 * (28. * Y00 - 42. * Y01 + 46. * Y02 + 6. * Y11 -
                         42. * Y12 + 28. * Y22)) +
             3. * (-2. * BB *
                       (46. * Y02 * Y11 - 14. * Y0 * Y111 + Y1111 -
                        14. * Y111 * Y2 - 14. * Y01 * Y22 + 28. * Y11 * Y22 +
                        Y00 * (28. * Y11 - 14. * Y12 + Y22)) +
                   5. * Z00 + 40. * Z01p - 2. * Z02 + 8. * Z1p1p + 40. * Z1p2 +
                   5. * Z22),
         -2. * AA *
                 (3. * X02 * X11 - 7. * X0 * X111 + 3. * X1111 -
                  7. * X111 * X2) -
             6. * BB * Y02 * Y11 + 14. * BB * Y0 * Y111 - 6. * BB * Y1111 -
             2. * AB *
                 (2. * X12 * Y01 - 7. * X01 * Y11 + X02 * Y11 - 7. * X12 * Y11 +
                  X11 * (-7. * Y01 + Y02 + 6. * Y11 - 7. * Y12) +
                  2. * X01 * Y12) +
             14. * BB * Y111 * Y2 - 5. * Z01p + Z02 - 7. * Z1p1p - 5. * Z1p2,
         -(AA * X1111) - 2. * AB * X11 * Y11 - BB * Y1111 + Z1p1p});
    // Update Volume
    auto m0 = scalar_type(0.0);
    for (size_t i = 0; i < 3; ++i) {
      m0 += coeffs[i] * m0_basis[i];
    }
    // Update Centroid
    auto m1x = scalar_type(0.0), m1y = scalar_type(0.0), m1z = scalar_type(0.0);
    for (size_t i = 0; i < 4; ++i) {
      m1x += coeffs[3 + i] * m1x_basis[i];
      m1y += coeffs[3 + i] * m1y_basis[i];
    }
    for (size_t i = 0; i < 5; ++i) {
      m1z += coeffs[7 + i] * m1z_basis[i];
    }
    // m0 *= area_proj_triangle;
    m1x *= area_proj_triangle;
    m1y *= area_proj_triangle;
    m1z *= area_proj_triangle;
    moments_with_gradient.volume() = m0.value();
    moments_with_gradient.volume_gradient() = m0.gradient();
    Pt& centroid = moments_with_gradient.centroid().getPt();
    auto& centroid_gradient = moments_with_gradient.centroid().getData();
    centroid[0] = m1x.value();
    centroid[1] = m1y.value();
    centroid[2] = m1z.value();
    centroid_gradient[0] = m1x.gradient();
    centroid_gradient[1] = m1y.gradient();
    centroid_gradient[2] = m1z.gradient();
  }

  return moments_with_gradient;
}  // namespace IRL

template <>
inline enable_if_t<!has_embedded_gradient<VolumeMoments>::value, VolumeMoments>
computeType3Contribution(const AlignedParaboloid& a_paraboloid,
                         const RationalBezierArc& a_arc) {
  auto moments = VolumeMoments::fromScalarConstant(0.0);
  const auto& pt_0 = a_arc.start_point();
  const auto& cp = a_arc.control_point();
  const auto& pt_1 = a_arc.end_point();
  const auto weight = a_arc.weight();
  const auto A = a_paraboloid.a(), B = a_paraboloid.b();
  const auto X0 = pt_0[0], X1 = cp[0], X2 = pt_1[0];
  const auto Y0 = pt_0[1], Y1 = cp[1], Y2 = pt_1[1];
  const auto Z0 = pt_0[2], Z1p = cp[2], Z2 = pt_1[2];
  const double Z1P = -A * X1 * X1 - B * Y1 * Y1;
  const double AA = A * A, BB = B * B, AB = A * B;
  const double X00 = X0 * X0, X11 = X1 * X1, X22 = X2 * X2;
  const double X000 = X00 * X0, X111 = X11 * X1, X222 = X22 * X2;
  const double X0000 = X00 * X00, X1111 = X11 * X11, X2222 = X22 * X22;
  const double Y00 = Y0 * Y0, Y11 = Y1 * Y1, Y22 = Y2 * Y2;
  const double Y000 = Y00 * Y0, Y111 = Y11 * Y1, Y222 = Y22 * Y2;
  const double Y0000 = Y00 * Y00, Y1111 = Y11 * Y11, Y2222 = Y22 * Y22;
  const double Z00 = Z0 * Z0, Z22 = Z2 * Z2;
  const double Z1p1p = Z1p * Z1p, Z1P1P = Z1P * Z1P;
  const double X02 = X0 * X2, X12 = X1 * X2, X01 = X0 * X1;
  const double Y02 = Y0 * Y2, Y12 = Y1 * Y2, Y01 = Y0 * Y1;
  const double Z02 = Z0 * Z2, Z01p = Z0 * Z1p, Z01P = Z0 * Z1P;
  const double Z1p1P = Z1p * Z1P, Z1p2 = Z1p * Z2, Z1P2 = Z1P * Z2;
  const double X0Z0 = X0 * Z0, X0Z1p = X0 * Z1p, X0Z1P = X0 * Z1P,
               X0Z2 = X0 * Z2;
  const double X1Z0 = X1 * Z0, X1Z1p = X1 * Z1p, X1Z1P = X1 * Z1P,
               X1Z2 = X1 * Z2;
  const double X2Z0 = X2 * Z0, X2Z1p = X2 * Z1p, X2Z1P = X2 * Z1P,
               X2Z2 = X2 * Z2;
  const double Y0Z0 = Y0 * Z0, Y0Z1p = Y0 * Z1p, Y0Z1P = Y0 * Z1P,
               Y0Z2 = Y0 * Z2;
  const double Y1Z0 = Y1 * Z0, Y1Z1p = Y1 * Z1p, Y1Z1P = Y1 * Z1P,
               Y1Z2 = Y1 * Z2;
  const double Y2Z0 = Y2 * Z0, Y2Z1p = Y2 * Z1p, Y2Z1P = Y2 * Z1P,
               Y2Z2 = Y2 * Z2;
  const double area_proj_triangle =
      0.5 * (X0 * (Y2 - Y1) + X2 * (Y1 - Y0) + X1 * (Y0 - Y2));
  assert(weight >= 0.0);
  // Compute coefficients (functions of the weight)
  std::array<double, 12> coeffs;
  if (weight < 0.35)  // We use the exact expressions
  {
    coeffs = coeffsV3andC3Exact<double>(weight);
  } else if (weight < 1.7)  // We use the 40th order Taylor series (w -> 1)
  {
    coeffs = coeffsV3andC3SeriesOne<double>(weight);
  } else if (weight < 1.0e9)  // We use the exact expressions
  {
    coeffs = coeffsV3andC3Exact<double>(weight);
  }
  // else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
  //   coeffs = coeffsV3andC3SeriesInfinity(weight);
  else  // This is within DBL_EPSILON of the actual value
  {
    coeffs = std::array<double, 12>({1.0 / 6.0, 2.0 / 3.0, 0.0, -1.0 / 140.0,
                                     1.0 / 420.0, -1.0 / 280.0, 1.0 / 14.0,
                                     1.0 / 630.0, -1.0 / 1890.0, 1.0 / 3780.0,
                                     -1.0 / 252.0, 1.0 / 18.0});
  }
  // Compute basis (combinations of pt_0, cp, and pt_1)
  auto m0_basis = std::array<double, 3>(
      {signedDistance(0.5 * (pt_0 + pt_1), a_paraboloid),
       signedDistance(0.25 * (pt_0 + pt_1) + 0.5 * cp, a_paraboloid),
       signedDistance(cp, a_paraboloid)});
  auto m1x_basis = std::array<double, 4>(
      {-6. * (X0Z0 - X0Z2 - X2Z0 + X2Z2 - 2. * B * X2 * Y00 +
              2. * B * X0 * Y02 + 2. * B * X2 * Y02 - 2. * B * X0 * Y22),
       2. * (5. * X0Z0 + 10. * X0Z1p + 6. * X0Z1P + 7. * X0Z2 +
             30. * A * X02 * X1 - 11. * X1Z0 - 4. * X1Z1p - 11. * X1Z2 +
             7. * X2Z0 + 10. * X2Z1p + 6. * X2Z1P + 5. * X2Z2 -
             14. * B * X1 * Y00 + 4. * B * X2 * Y00 + 14. * B * X0 * Y01 -
             4. * B * X1 * Y01 + 10. * B * X2 * Y01 - 4. * B * X0 * Y02 +
             10. * B * X1 * Y02 - 4. * B * X2 * Y02 + 4. * B * X0 * Y11 +
             4. * B * X2 * Y11 + 10. * B * X0 * Y12 - 4. * B * X1 * Y12 +
             14. * B * X2 * Y12 + 4. * B * X0 * Y22 - 14. * B * X1 * Y22),
       2. * (-5. * X0Z1p + 18. * X0Z1P + X0Z2 + 6. * A * X02 * X1 - 5. * X1Z0 -
             6. * X1Z1p - 6. * X1Z1P - 5. * X1Z2 + X2Z0 - 5. * X2Z1p +
             18. * X2Z1P - 12. * B * X1 * Y01 + 2. * B * X2 * Y01 +
             2. * B * X1 * Y02 + 12. * B * X0 * Y11 + 12. * B * X2 * Y11 +
             2. * B * X0 * Y12 - 12. * B * X1 * Y12),
       2. * (X1Z1p - X1Z1P)});
  auto m1y_basis = std::array<double, 4>(
      {6. * (-Y0Z0 + Y0Z2 + 2. * A * (X22 * Y0 + X00 * Y2 - X02 * (Y0 + Y2)) +
             Y2Z0 - Y2Z2),
       2. *
           (5. * Y0Z0 + 10. * Y0Z1p + 6. * Y0Z1P + 7. * Y0Z2 +
            30. * B * Y02 * Y1 - 11. * Y1Z0 - 4. * Y1Z1p - 11. * Y1Z2 +
            2. * A *
                (-2. * X02 * Y0 + 2. * X11 * Y0 + 5. * X12 * Y0 +
                 2. * X22 * Y0 - 7. * X00 * Y1 + 5. * X02 * Y1 - 2. * X12 * Y1 -
                 7. * X22 * Y1 + 2. * X00 * Y2 - 2. * X02 * Y2 + 2. * X11 * Y2 +
                 7. * X12 * Y2 + X01 * (7. * Y0 - 2. * Y1 + 5. * Y2)) +
            7. * Y2Z0 + 10. * Y2Z1p + 6. * Y2Z1P + 5. * Y2Z2),
       -2. * (5. * Y0Z1p - 18. * Y0Z1P - Y0Z2 - 6. * B * Y02 * Y1 + 5. * Y1Z0 +
              6. * Y1Z1p + 6. * Y1Z1P + 5. * Y1Z2 -
              2. * A *
                  (X12 * Y0 - 6. * X01 * Y1 + X02 * Y1 - 6. * X12 * Y1 +
                   X01 * Y2 + 6. * X11 * (Y0 + Y2)) -
              Y2Z0 + 5. * Y2Z1p - 18. * Y2Z1P),
       2. * (Y1Z1p - Y1Z1P)});
  auto m1z_basis = std::array<double, 5>(
      {-(AA * (21. * X0000 + 28. * X000 * X2 + 30. * X00 * X22 +
               28. * X0 * X222 + 21. * X2222)) -
           21. * BB * Y0000 - 28. * BB * Y000 * Y2 - 30. * BB * Y00 * Y22 -
           2. * AB *
               (X00 * (21. * Y00 + 14. * Y02 + 5. * Y22) +
                2. * X02 * (7. * Y00 + 10. * Y02 + 7. * Y22) +
                X22 * (5. * Y00 + 14. * Y02 + 21. * Y22)) -
           28. * BB * Y0 * Y222 - 21. * BB * Y2222 + 40. * Z00 + 48. * Z02 +
           40. * Z22,
       3. * AA *
               (21. * X000 * X1 - 7. * X00 * X11 - 10. * X02 * X11 +
                35. * X00 * X12 - 7. * X000 * X2 - 10. * X00 * X22 +
                35. * X01 * X22 - 7. * X11 * X22 - 7. * X0 * X222 +
                21. * X1 * X222) -
           AB * (7. * X11 * Y00 - 35. * X12 * Y00 + 10. * X22 * Y00 -
                 63. * X00 * Y01 + 20. * X12 * Y01 - 35. * X22 * Y01 +
                 21. * X00 * Y02 + 10. * X11 * Y02 - 70. * X12 * Y02 +
                 21. * X22 * Y02 + 7. * X00 * Y11 + 7. * X22 * Y11 -
                 35. * X00 * Y12 + 28. * X12 * Y12 - 63. * X22 * Y12 +
                 X01 * (-63. * Y00 + 28. * Y01 - 70. * Y02 + 20. * Y12 -
                        35. * Y22) +
                 10. * X00 * Y22 + 7. * X11 * Y22 - 63. * X12 * Y22 +
                 X02 * (21. * Y00 - 70. * Y01 + 40. * Y02 + 10. * Y11 -
                        70. * Y12 + 21. * Y22)) -
           3. * (BB * (7. * Y00 * Y11 + 10. * Y02 * Y11 - 35. * Y00 * Y12 +
                       7. * Y000 * (-3. * Y1 + Y2) + 10. * Y00 * Y22 -
                       35. * Y01 * Y22 + 7. * Y11 * Y22 + 7. * Y0 * Y222 -
                       21. * Y1 * Y222) +
                 2. * (5. * Z00 + 10. * Z01p + 4. * Z02 - 2. * Z1p1p +
                       10. * Z1p2 + 5. * Z22)),
       -6. * AA *
               (46. * X02 * X11 - 14. * X0 * X111 + X1111 - 14. * X111 * X2 -
                14. * X01 * X22 + 28. * X11 * X22 +
                X00 * (28. * X11 - 14. * X12 + X22)) -
           2. * AB *
               (X22 * Y00 + 112. * X01 * Y01 - 28. * X02 * Y01 -
                14. * X22 * Y01 - 28. * X01 * Y02 + 4. * X02 * Y02 +
                28. * X00 * Y11 - 42. * X01 * Y11 + 46. * X02 * Y11 +
                28. * X22 * Y11 -
                2. * X12 *
                    (7. * Y00 - 46. * Y01 + 14. * Y02 + 21. * Y11 - 56. * Y12) -
                14. * X00 * Y12 + 92. * X01 * Y12 - 28. * X02 * Y12 +
                X00 * Y22 - 14. * X01 * Y22 +
                X11 * (28. * Y00 - 42. * Y01 + 46. * Y02 + 6. * Y11 -
                       42. * Y12 + 28. * Y22)) +
           3. * (-2. * BB *
                     (46. * Y02 * Y11 - 14. * Y0 * Y111 + Y1111 -
                      14. * Y111 * Y2 - 14. * Y01 * Y22 + 28. * Y11 * Y22 +
                      Y00 * (28. * Y11 - 14. * Y12 + Y22)) +
                 5. * Z00 + 40. * Z01p - 2. * Z02 + 8. * Z1p1p + 40. * Z1p2 +
                 5. * Z22),
       -2. * AA *
               (3. * X02 * X11 - 7. * X0 * X111 + 3. * X1111 - 7. * X111 * X2) -
           6. * BB * Y02 * Y11 + 14. * BB * Y0 * Y111 - 6. * BB * Y1111 -
           2. * AB *
               (2. * X12 * Y01 - 7. * X01 * Y11 + X02 * Y11 - 7. * X12 * Y11 +
                X11 * (-7. * Y01 + Y02 + 6. * Y11 - 7. * Y12) +
                2. * X01 * Y12) +
           14. * BB * Y111 * Y2 - 5. * Z01p + Z02 - 7. * Z1p1p - 5. * Z1p2,
       -(AA * X1111) - 2. * AB * X11 * Y11 - BB * Y1111 + Z1p1p});
  for (size_t i = 0; i < 3; ++i) {
    moments.volume() += coeffs[i] * m0_basis[i];
  }
  for (size_t i = 0; i < 4; ++i) {
    moments.centroid()[0] += coeffs[3 + i] * m1x_basis[i];
    moments.centroid()[1] += coeffs[3 + i] * m1y_basis[i];
  }
  for (size_t i = 0; i < 5; ++i) {
    moments.centroid()[2] += coeffs[7 + i] * m1z_basis[i];
  }
  moments.volume() *= area_proj_triangle;
  moments.centroid()[0] *= area_proj_triangle;
  moments.centroid()[1] *= area_proj_triangle;
  moments.centroid()[2] *= area_proj_triangle;
  // std::cout << "M3 contribution = " << moments.volume()
  //           << " ------- ARC = " << pt_0 << " -- " << cp << " -- " << pt_1
  //           << " -- w = " << weight << std::endl;
  return moments;
}  // namespace IRL

template <class ReturnType, class PtType>
ReturnType computeFaceOnlyContribution(const AlignedParaboloid& a_paraboloid,
                                       const Plane& a_face_plane,
                                       const PtType& a_pt_ref);

template <>
inline enable_if_t<!has_embedded_gradient<Volume>::value, Volume>
computeFaceOnlyContribution(const AlignedParaboloid& a_paraboloid,
                            const Plane& a_face_plane, const Pt& a_pt_ref) {
  assert(a_paraboloid.a() * a_paraboloid.b() > 0.0);
  assert(std::fabs(a_face_plane.normal()[2]) > DBL_EPSILON);
  const double a =
      -a_face_plane.normal()[0] / safelyEpsilon(a_face_plane.normal()[2]);
  const double b =
      -a_face_plane.normal()[1] / safelyEpsilon(a_face_plane.normal()[2]);
  const double c =
      a_face_plane.distance() / safelyEpsilon(a_face_plane.normal()[2]);
  const double factor = 4.0 * a_paraboloid.a() * a_paraboloid.b() * c -
                        a_paraboloid.a() * b * b - a_paraboloid.b() * a * a;
  return std::copysign(
      M_PI * factor * factor /
          (32.0 * std::pow(a_paraboloid.a() * a_paraboloid.b(), 2.5)),
      -a_face_plane.normal()[2]);
}

template <>
inline enable_if_t<!has_embedded_gradient<VolumeMoments>::value, VolumeMoments>
computeFaceOnlyContribution(const AlignedParaboloid& a_paraboloid,
                            const Plane& a_face_plane, const Pt& a_pt_ref) {
  assert(a_paraboloid.a() * a_paraboloid.b() > 0.0);
  assert(std::fabs(a_face_plane.normal()[2]) > DBL_EPSILON);
  const double a = -a_face_plane.normal()[0] / a_face_plane.normal()[2];
  const double b = -a_face_plane.normal()[1] / a_face_plane.normal()[2];
  const double c = a_face_plane.distance() / a_face_plane.normal()[2];
  const auto A = a_paraboloid.a(), B = a_paraboloid.b();
  auto moments = VolumeMoments::fromScalarConstant(0.0);
  const double factor = (a * a * B + A * (b * b - 4. * B * c)) *
                        (a * a * B + A * (b * b - 4. * B * c)) * M_PI;
  moments.volume() = std::copysign(factor / (32. * std::pow(A * B, 2.5)),
                                   -a_face_plane.normal()[2]);
  moments.centroid()[0] = a * B *
                          std::copysign(factor / (64. * std::pow(A * B, 3.5)),
                                        a_face_plane.normal()[2]);
  moments.centroid()[1] = b * A *
                          std::copysign(factor / (64. * std::pow(A * B, 3.5)),
                                        a_face_plane.normal()[2]);
  moments.centroid()[2] =
      (5. * A * (b * b) + 5. * (a * a) * B - 8. * A * B * c) *
      std::copysign(factor / (384. * std::pow(A * B, 3.5)),
                    a_face_plane.normal()[2]);
  return moments;
}

template <class ReturnType>
inline enable_if_t<has_embedded_gradient<ReturnType>::value, ReturnType>
computeFaceOnlyContribution(
    const AlignedParaboloid& a_paraboloid, const Plane& a_face_plane,
    const PtWithGradient<typename ReturnType::gradient_type>& a_pt_ref) {
  using gradient_type = typename ReturnType::gradient_type;
  using scalar_type = ScalarWithGradient<gradient_type>;
  auto moments_with_gradient = ReturnType::fromScalarConstant(0.0);
  assert(a_paraboloid.a() * a_paraboloid.b() > 0.0);
  assert(std::fabs(a_face_plane.normal()[2]) > DBL_EPSILON);

  // Compute plane gradient
  const auto& pt_ref = a_pt_ref.getPt();
  const auto& pt_ref_grad = a_pt_ref.getData();
  auto face_normal = a_face_plane.normal();
  auto face_normal_withgrad = PtWithGradient<gradient_type>(
      Pt(face_normal[0], face_normal[1], face_normal[2]));
  auto& face_normal_grad = face_normal_withgrad.getData();
  face_normal_grad[0].setGradRx(0.0);
  face_normal_grad[1].setGradRx(face_normal[2]);
  face_normal_grad[2].setGradRx(-face_normal[1]);
  face_normal_grad[0].setGradRy(-face_normal[2]);
  face_normal_grad[1].setGradRy(0.0);
  face_normal_grad[2].setGradRy(face_normal[0]);
  face_normal_grad[0].setGradRz(face_normal[1]);
  face_normal_grad[1].setGradRz(-face_normal[0]);
  face_normal_grad[2].setGradRz(0.0);
  auto face_distance = face_normal[0] * pt_ref[0] + face_normal[1] * pt_ref[1] +
                       face_normal[2] * pt_ref[2];
  auto face_distance_grad =
      face_normal_grad[0] * pt_ref[0] + face_normal_grad[1] * pt_ref[1] +
      face_normal_grad[2] * pt_ref[2] + face_normal[0] * pt_ref_grad[0] +
      face_normal[1] * pt_ref_grad[1] + face_normal[2] * pt_ref_grad[2];

  const double a_value = -face_normal[0] / safelyEpsilon(face_normal[2]);
  const double b_value = -face_normal[1] / safelyEpsilon(face_normal[2]);
  const double c_value = face_distance / safelyEpsilon(face_normal[2]);
  const auto a_grad = -face_normal_grad[0] / safelyEpsilon(face_normal[2]) +
                      face_normal[0] * face_normal_grad[2] /
                          safelyEpsilon(face_normal[2] * face_normal[2]);
  const auto b_grad = -face_normal_grad[1] / safelyEpsilon(face_normal[2]) +
                      face_normal[1] * face_normal_grad[2] /
                          safelyEpsilon(face_normal[2] * face_normal[2]);
  const auto c_grad = face_distance_grad / safelyEpsilon(face_normal[2]) -
                      face_distance * face_normal_grad[2] /
                          safelyEpsilon(face_normal[2] * face_normal[2]);
  const auto a = scalar_type(a_value, a_grad);
  const auto b = scalar_type(b_value, b_grad);
  const auto c = scalar_type(c_value, c_grad);
  auto A_grad = gradient_type(0.0), B_grad = gradient_type(0.0);
  A_grad.setGradA(1.0);
  B_grad.setGradB(1.0);
  const auto A = scalar_type(a_paraboloid.a(), A_grad),
             B = scalar_type(a_paraboloid.b(), B_grad);
  const auto factor = (a * a * B + A * (b * b - 4. * B * c)) *
                      (a * a * B + A * (b * b - 4. * B * c)) * M_PI;
  const double sign_normal = std::copysign(1.0, a_face_plane.normal()[2]);
  const auto m0 =
      -sign_normal * factor / (32. * PowMoments<scalar_type>(A * B, 2.5));
  moments_with_gradient.volume() = m0.value();
  moments_with_gradient.volume_gradient() = m0.gradient();
  if constexpr (std::is_same<typename ReturnType::moment_type,
                             VolumeMoments>::value) {
    const auto m1x = sign_normal * a * B * factor /
                     (64. * PowMoments<scalar_type>(A * B, 3.5));
    const auto m1y = sign_normal * b * A * factor /
                     (64. * PowMoments<scalar_type>(A * B, 3.5));
    const auto m1z = sign_normal *
                     (5. * A * (b * b) + 5. * (a * a) * B - 8. * A * B * c) *
                     factor / (384. * PowMoments<scalar_type>(A * B, 3.5));
    Pt& centroid = moments_with_gradient.centroid().getPt();
    auto& centroid_gradient = moments_with_gradient.centroid().getData();
    centroid[0] = m1x.value();
    centroid[1] = m1y.value();
    centroid[2] = m1z.value();
    centroid_gradient[0] = m1x.gradient();
    centroid_gradient[1] = m1y.gradient();
    centroid_gradient[2] = m1z.gradient();
  }
  return moments_with_gradient;
}

template <class ReturnType, class PtType>
ReturnType computeTriangleCorrection(const AlignedParaboloid& a_paraboloid,
                                     const PtType& a_pt_0, const PtType& a_pt_1,
                                     const PtType& a_pt_2);

template <>
inline enable_if_t<!has_embedded_gradient<Volume>::value, Volume>
computeTriangleCorrection(const AlignedParaboloid& a_paraboloid,
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

template <class ReturnType>
inline enable_if_t<has_embedded_gradient<ReturnType>::value, ReturnType>
computeTriangleCorrection(
    const AlignedParaboloid& a_paraboloid,
    const PtWithGradient<typename ReturnType::gradient_type>& a_pt_0,
    const PtWithGradient<typename ReturnType::gradient_type>& a_pt_1,
    const PtWithGradient<typename ReturnType::gradient_type>& a_pt_2) {
  using gradient_type = typename ReturnType::gradient_type;
  using scalar_type = ScalarWithGradient<gradient_type>;
  auto moments_with_gradient = ReturnType::fromScalarConstant(0.0);
  std::array<scalar_type, 3> pt_0, pt_1, pt_2;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_0[d] = scalar_type(a_pt_0.getPt()[d], a_pt_0.getData()[d]);
    pt_1[d] = scalar_type(a_pt_1.getPt()[d], a_pt_1.getData()[d]);
    pt_2[d] = scalar_type(a_pt_2.getPt()[d], a_pt_2.getData()[d]);
  }
  auto A_grad = gradient_type(0.0), B_grad = gradient_type(0.0);
  A_grad.setGradA(1.0);
  B_grad.setGradB(1.0);
  const auto A = scalar_type(a_paraboloid.a(), A_grad),
             B = scalar_type(a_paraboloid.b(), B_grad);
  const auto X0 = pt_0[0], X1 = pt_1[0], X2 = pt_2[0];
  const auto Y0 = pt_0[1], Y1 = pt_1[1], Y2 = pt_2[1];
  const auto Z0 = pt_0[2], Z1 = pt_1[2], Z2 = pt_2[2];
  const auto triangle_area =
      0.5 * ((pt_1[1] - pt_2[1]) * pt_0[0] + (pt_2[1] - pt_0[1]) * pt_1[0] +
             (pt_0[1] - pt_1[1]) * pt_2[0]);
  const auto m0 = (-A * (pt_0[0] + pt_1[0]) * (pt_1[0] + pt_2[0]) +
                   -B * (pt_0[1] + pt_1[1]) * (pt_1[1] + pt_2[1]) - pt_0[2] -
                   2.0 * pt_1[2] - pt_2[2]) *
                  triangle_area / 6.0;
  moments_with_gradient.volume() = m0.value();
  moments_with_gradient.volume_gradient() = m0.gradient();
  if constexpr (std::is_same<typename ReturnType::moment_type,
                             VolumeMoments>::value) {
    const auto m1x =
        triangle_area *
        ((A * (-(X0 * X0 * X0) - X1 * X1 * X1 - X1 * X1 * X2 - X1 * (X2 * X2) -
               X2 * X2 * X2 - X0 * X0 * (X1 + X2) -
               X0 * (X1 * X1 + X1 * X2 + X2 * X2))) /
             10. +
         (B * (-(X1 * (Y0 * Y0 + 2. * Y0 * Y1 + 3. * (Y1 * Y1) + Y0 * Y2 +
                       2. * Y1 * Y2 + Y2 * Y2)) -
               X2 * (Y0 * Y0 + Y0 * Y1 + Y1 * Y1 + 2. * Y0 * Y2 + 2. * Y1 * Y2 +
                     3. * (Y2 * Y2)) -
               X0 * (3. * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
                     2. * Y0 * (Y1 + Y2)))) /
             30. +
         (-(X0 * (2. * Z0 + Z1 + Z2)) - X1 * (Z0 + 2. * Z1 + Z2) -
          X2 * (Z0 + Z1 + 2. * Z2)) /
             12.);
    const auto m1y =
        -triangle_area *
        ((B * (Y0 * Y0 * Y0 + Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) +
               Y2 * Y2 * Y2 + Y0 * Y0 * (Y1 + Y2) +
               Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2))) /
             10. +
         (A * (X0 * X0 * (3. * Y0 + Y1 + Y2) + X1 * X1 * (Y0 + 3. * Y1 + Y2) +
               X2 * X2 * (Y0 + Y1 + 3. * Y2) + X1 * X2 * (Y0 + 2. * (Y1 + Y2)) +
               X0 * (X1 * (2. * Y0 + 2. * Y1 + Y2) +
                     X2 * (2. * Y0 + Y1 + 2. * Y2)))) /
             30. +
         (Y0 * (2. * Z0 + Z1 + Z2) + Y1 * (Z0 + 2. * Z1 + Z2) +
          Y2 * (Z0 + Z1 + 2. * Z2)) /
             12.);
    const auto m1z =
        triangle_area *
        ((A * A *
          (X0 * X0 * X0 * X0 + X1 * X1 * X1 * X1 + X1 * X1 * X1 * X2 +
           X1 * X1 * (X2 * X2) + X1 * (X2 * X2 * X2) + X2 * X2 * X2 * X2 +
           X0 * X0 * X0 * (X1 + X2) + X0 * X0 * (X1 * X1 + X1 * X2 + X2 * X2) +
           X0 *
               (X1 * X1 * X1 + X1 * X1 * X2 + X1 * (X2 * X2) + X2 * X2 * X2))) /
             30. +
         (B * B *
          (Y0 * Y0 * Y0 * Y0 + Y1 * Y1 * Y1 * Y1 + Y1 * Y1 * Y1 * Y2 +
           Y1 * Y1 * (Y2 * Y2) + Y1 * (Y2 * Y2 * Y2) + Y2 * Y2 * Y2 * Y2 +
           Y0 * Y0 * Y0 * (Y1 + Y2) + Y0 * Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2) +
           Y0 *
               (Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) + Y2 * Y2 * Y2))) /
             30. +
         (A * B *
          (X1 * X2 *
               (Y0 * Y0 + 3. * (Y1 * Y1) + 4. * Y1 * Y2 + 3. * (Y2 * Y2) +
                2. * Y0 * (Y1 + Y2)) +
           X0 * X0 *
               (6. * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
                3. * Y0 * (Y1 + Y2)) +
           X1 * X1 *
               (Y0 * Y0 + 6. * (Y1 * Y1) + 3. * Y1 * Y2 + Y2 * Y2 +
                Y0 * (3. * Y1 + Y2)) +
           X2 * X2 *
               (Y0 * Y0 + Y1 * Y1 + 3. * Y1 * Y2 + 6. * (Y2 * Y2) +
                Y0 * (Y1 + 3. * Y2)) +
           X0 * (X1 * (3. * (Y0 * Y0) + 4. * Y0 * Y1 + 3. * (Y1 * Y1) +
                       2. * Y0 * Y2 + 2. * Y1 * Y2 + Y2 * Y2) +
                 X2 * (3. * (Y0 * Y0) + 2. * Y0 * Y1 + Y1 * Y1 + 4. * Y0 * Y2 +
                       2. * Y1 * Y2 + 3. * (Y2 * Y2))))) /
             90. +
         (-(Z0 * Z0) - Z1 * Z1 - Z1 * Z2 - Z2 * Z2 - Z0 * (Z1 + Z2)) / 12.);
    Pt& centroid = moments_with_gradient.centroid().getPt();
    auto& centroid_gradient = moments_with_gradient.centroid().getData();
    centroid[0] = m1x.value();
    centroid[1] = m1y.value();
    centroid[2] = m1z.value();
    centroid_gradient[0] = m1x.gradient();
    centroid_gradient[1] = m1y.gradient();
    centroid_gradient[2] = m1z.gradient();
  }
  return moments_with_gradient;
}

template <>
inline enable_if_t<!has_embedded_gradient<VolumeMoments>::value, VolumeMoments>
computeTriangleCorrection(const AlignedParaboloid& a_paraboloid,
                          const Pt& a_pt_0, const Pt& a_pt_1,
                          const Pt& a_pt_2) {
  auto moments = VolumeMoments::fromScalarConstant(0.0);
  const double A = a_paraboloid.a(), B = a_paraboloid.b();
  const double X0 = a_pt_0[0], X1 = a_pt_1[0], X2 = a_pt_2[0];
  const double Y0 = a_pt_0[1], Y1 = a_pt_1[1], Y2 = a_pt_2[1];
  const double Z0 = a_pt_0[2], Z1 = a_pt_1[2], Z2 = a_pt_2[2];
  const double triangle_area = 0.5 * ((a_pt_1[1] - a_pt_2[1]) * a_pt_0[0] +
                                      (a_pt_2[1] - a_pt_0[1]) * a_pt_1[0] +
                                      (a_pt_0[1] - a_pt_1[1]) * a_pt_2[0]);
  moments.volume() =
      (-a_paraboloid.a() * (a_pt_0[0] + a_pt_1[0]) * (a_pt_1[0] + a_pt_2[0]) +
       -a_paraboloid.b() * (a_pt_0[1] + a_pt_1[1]) * (a_pt_1[1] + a_pt_2[1]) -
       a_pt_0[2] - 2.0 * a_pt_1[2] - a_pt_2[2]) *
      triangle_area / 6.0;
  moments.centroid()[0] =
      triangle_area *
      ((A * (-(X0 * X0 * X0) - X1 * X1 * X1 - X1 * X1 * X2 - X1 * (X2 * X2) -
             X2 * X2 * X2 - X0 * X0 * (X1 + X2) -
             X0 * (X1 * X1 + X1 * X2 + X2 * X2))) /
           10. +
       (B * (-(X1 * (Y0 * Y0 + 2. * Y0 * Y1 + 3. * (Y1 * Y1) + Y0 * Y2 +
                     2. * Y1 * Y2 + Y2 * Y2)) -
             X2 * (Y0 * Y0 + Y0 * Y1 + Y1 * Y1 + 2. * Y0 * Y2 + 2. * Y1 * Y2 +
                   3. * (Y2 * Y2)) -
             X0 * (3. * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
                   2. * Y0 * (Y1 + Y2)))) /
           30. +
       (-(X0 * (2. * Z0 + Z1 + Z2)) - X1 * (Z0 + 2. * Z1 + Z2) -
        X2 * (Z0 + Z1 + 2. * Z2)) /
           12.);
  moments.centroid()[1] =
      -triangle_area *
      ((B * (Y0 * Y0 * Y0 + Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) +
             Y2 * Y2 * Y2 + Y0 * Y0 * (Y1 + Y2) +
             Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2))) /
           10. +
       (A * (X0 * X0 * (3. * Y0 + Y1 + Y2) + X1 * X1 * (Y0 + 3. * Y1 + Y2) +
             X2 * X2 * (Y0 + Y1 + 3. * Y2) + X1 * X2 * (Y0 + 2. * (Y1 + Y2)) +
             X0 * (X1 * (2. * Y0 + 2. * Y1 + Y2) +
                   X2 * (2. * Y0 + Y1 + 2. * Y2)))) /
           30. +
       (Y0 * (2. * Z0 + Z1 + Z2) + Y1 * (Z0 + 2. * Z1 + Z2) +
        Y2 * (Z0 + Z1 + 2. * Z2)) /
           12.);
  moments.centroid()[2] =
      triangle_area *
      ((A * A *
        (X0 * X0 * X0 * X0 + X1 * X1 * X1 * X1 + X1 * X1 * X1 * X2 +
         X1 * X1 * (X2 * X2) + X1 * (X2 * X2 * X2) + X2 * X2 * X2 * X2 +
         X0 * X0 * X0 * (X1 + X2) + X0 * X0 * (X1 * X1 + X1 * X2 + X2 * X2) +
         X0 * (X1 * X1 * X1 + X1 * X1 * X2 + X1 * (X2 * X2) + X2 * X2 * X2))) /
           30. +
       (B * B *
        (Y0 * Y0 * Y0 * Y0 + Y1 * Y1 * Y1 * Y1 + Y1 * Y1 * Y1 * Y2 +
         Y1 * Y1 * (Y2 * Y2) + Y1 * (Y2 * Y2 * Y2) + Y2 * Y2 * Y2 * Y2 +
         Y0 * Y0 * Y0 * (Y1 + Y2) + Y0 * Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2) +
         Y0 * (Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) + Y2 * Y2 * Y2))) /
           30. +
       (A * B *
        (X1 * X2 *
             (Y0 * Y0 + 3. * (Y1 * Y1) + 4. * Y1 * Y2 + 3. * (Y2 * Y2) +
              2. * Y0 * (Y1 + Y2)) +
         X0 * X0 *
             (6. * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
              3. * Y0 * (Y1 + Y2)) +
         X1 * X1 *
             (Y0 * Y0 + 6. * (Y1 * Y1) + 3. * Y1 * Y2 + Y2 * Y2 +
              Y0 * (3. * Y1 + Y2)) +
         X2 * X2 *
             (Y0 * Y0 + Y1 * Y1 + 3. * Y1 * Y2 + 6. * (Y2 * Y2) +
              Y0 * (Y1 + 3. * Y2)) +
         X0 * (X1 * (3. * (Y0 * Y0) + 4. * Y0 * Y1 + 3. * (Y1 * Y1) +
                     2. * Y0 * Y2 + 2. * Y1 * Y2 + Y2 * Y2) +
               X2 * (3. * (Y0 * Y0) + 2. * Y0 * Y1 + Y1 * Y1 + 4. * Y0 * Y2 +
                     2. * Y1 * Y2 + 3. * (Y2 * Y2))))) /
           90. +
       (-(Z0 * Z0) - Z1 * Z1 - Z1 * Z2 - Z2 * Z2 - Z0 * (Z1 + Z2)) / 12.);
  return moments;
}

}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_MOMENT_CONTRIBUTIONS_TPP_
