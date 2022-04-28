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
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/helpers/mymath.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

namespace IRL {

inline double signedDistance(const Pt& a_pt,
                             const AlignedParaboloid& a_paraboloid) {
  return a_paraboloid.a() * a_pt[0] * a_pt[0] +
         a_paraboloid.b() * a_pt[1] * a_pt[1] + a_pt[2];
}

template <class ReturnType, class PtType>
ReturnType computeType1Contribution(const PtType& a_ref_pt,
                                    const PtType& a_pt_0, const PtType& a_pt_1);

template <>
inline enable_if_t<!has_embedded_gradient<Volume>::value, Volume>
computeType1Contribution(const Pt& a_ref_pt, const Pt& a_pt_0,
                         const Pt& a_pt_1) {
  return (a_ref_pt[2] + a_pt_0[2] + a_pt_1[2]) / 6.0 *
         ((a_pt_0[0] - a_ref_pt[0]) * (a_pt_1[1] - a_ref_pt[1]) -
          (a_pt_1[0] - a_ref_pt[0]) * (a_pt_0[1] - a_ref_pt[1]));
}

template <class ReturnType>
inline enable_if_t<has_embedded_gradient<ReturnType>::value, ReturnType>
computeType1Contribution(
    const PtWithGradient<typename ReturnType::gradient_type>& a_ref_pt,
    const PtWithGradient<typename ReturnType::gradient_type>& a_pt_0,
    const PtWithGradient<typename ReturnType::gradient_type>& a_pt_1) {
  auto volume_with_gradient = ReturnType::fromScalarConstant(0.0);
  const auto& pt_0 = a_pt_0.getPt();
  const auto& pt_1 = a_pt_1.getPt();
  const auto& ref_pt = a_ref_pt.getPt();
  const auto& pt_0_grad = a_pt_0.getData();
  const auto& pt_1_grad = a_pt_1.getData();
  const auto& ref_pt_grad = a_ref_pt.getData();
  volume_with_gradient.volume() =
      (ref_pt[2] + pt_0[2] + pt_1[2]) / 6.0 *
      ((pt_0[0] - ref_pt[0]) * (pt_1[1] - ref_pt[1]) -
       (pt_1[0] - ref_pt[0]) * (pt_0[1] - ref_pt[1]));
  volume_with_gradient.gradient() =
      (ref_pt_grad[2] + pt_0_grad[2] + pt_1_grad[2]) / 6.0 *
          ((pt_0[0] - ref_pt[0]) * (pt_1[1] - ref_pt[1]) -
           (pt_1[0] - ref_pt[0]) * (pt_0[1] - ref_pt[1])) +
      (ref_pt[2] + pt_0[2] + pt_1[2]) / 6.0 *
          ((pt_0_grad[0] - ref_pt_grad[0]) * (pt_1[1] - ref_pt[1]) +
           (pt_0[0] - ref_pt[0]) * (pt_1_grad[1] - ref_pt_grad[1]) -
           (pt_1_grad[0] - ref_pt_grad[0]) * (pt_0[1] - ref_pt[1]) -
           (pt_1[0] - ref_pt[0]) * (pt_0_grad[1] - ref_pt_grad[1]));
  return volume_with_gradient;
}

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
  auto volume_with_gradient = ReturnType::fromScalarConstant(0.0);
  const auto& pt_0 = a_pt_0.getPt();
  const auto& pt_1 = a_pt_1.getPt();
  const auto& pt_0_grad = a_pt_0.getData();
  const auto& pt_1_grad = a_pt_1.getData();
  const double A = a_aligned_paraboloid.a(), B = a_aligned_paraboloid.b();
  auto A_grad = gradient_type(0.0), B_grad = gradient_type(0.0);
  A_grad.setGradA(1.0);
  B_grad.setGradB(1.0);
  volume_with_gradient.volume() =
      (pt_0[0] * pt_1[1] - pt_1[0] * pt_0[1]) / 12.0 *
      (-pt_0[2] - pt_1[2] + A * pt_0[0] * pt_1[0] + B * pt_0[1] * pt_1[1]);
  volume_with_gradient.gradient() =
      (pt_0_grad[0] * pt_1[1] - pt_1[0] * pt_0_grad[1] +
       pt_0[0] * pt_1_grad[1] - pt_1_grad[0] * pt_0[1]) /
          12.0 *
          (-pt_0[2] - pt_1[2] + A * pt_0[0] * pt_1[0] + B * pt_0[1] * pt_1[1]) +
      (pt_0[0] * pt_1[1] - pt_1[0] * pt_0[1]) / 12.0 *
          (-pt_0_grad[2] - pt_1_grad[2] + A_grad * pt_0[0] * pt_1[0] +
           B_grad * pt_0[1] * pt_1[1] + A * pt_0_grad[0] * pt_1[0] +
           B * pt_0_grad[1] * pt_1[1] + A * pt_0[0] * pt_1_grad[0] +
           B * pt_0[1] * pt_1_grad[1]);
  return volume_with_gradient;
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

inline std::array<double, 3> coeffsV3SeriesOne(const double a_weight) {
  std::array<double, 3> coeffs;
  coeffs.fill(0.0);
  double x = 1.0;
  for (UnsignedIndex_t i = 0; i <= 40; ++i) {
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      coeffs[j] += v3Series[i][j] * x;
    }
    x *= a_weight - 1.0;
  }
  return coeffs;
}

template <class GradientType>
inline std::array<std::pair<double, GradientType>, 3>
coeffsV3SeriesOneWithGradient(const double a_weight,
                              const GradientType& a_weight_grad) {
  using dbl_gradient = typename std::pair<double, GradientType>;
  std::array<dbl_gradient, 3> coeffs;
  coeffs.fill(dbl_gradient({0.0, 0.0}));
  double x = 1.0;
  for (UnsignedIndex_t i = 0; i <= 40; ++i) {
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      coeffs[j].first += v3Series[i][j] * x;
    }
    x *= a_weight - 1.0;
  }
  x = 1.0;
  for (UnsignedIndex_t i = 1; i <= 40; ++i) {
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      coeffs[j].second +=
          static_cast<double>(i) * a_weight_grad * v3Series[i][j] * x;
    }
    x *= a_weight - 1.0;
  }
  return coeffs;
}

inline std::array<double, 12> coeffsV3andC3SeriesOne(const double a_weight) {
  std::array<double, 12> coeffs;
  coeffs.fill(0.0);
  double x = 1.0;
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

inline std::array<double, 12> coeffsV3andC3SeriesInfinity(
    const double a_weight) {
  std::cout << "Not implemented yet" << std::endl;
  std::exit(-1);
  return std::array<double, 12>();
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

template <class GradientType>
inline std::array<std::pair<double, GradientType>, 3> coeffsV3ExactWithGradient(
    const double a_weight, const GradientType& a_weight_grad) {
  using dbl_gradient = typename std::pair<double, GradientType>;
  const double w2 = a_weight * a_weight;
  const double w3 = w2 * a_weight;
  const double w4 = w2 * w2;
  const double w5 = w2 * w3;
  const double w6 = w3 * w3;
  const GradientType w2_grad = 2.0 * a_weight_grad * a_weight;
  const GradientType w3_grad = 3.0 * a_weight_grad * w2;
  const GradientType w4_grad = 4.0 * a_weight_grad * w3;
  const GradientType w5_grad = 5.0 * a_weight_grad * w4;
  const GradientType w6_grad = 6.0 * a_weight_grad * w5;
  const double L = 1.0 / ((a_weight - 1.0) * (a_weight + 1.0));
  const GradientType L_grad = -2.0 * a_weight_grad * a_weight * L * L;
  const double L3 = L * L * L;
  const GradientType L3_grad = 3.0 * L_grad * L * L;
  const double S = (a_weight < 1.0) ? sqrt(1.0 - a_weight * a_weight)
                                    : sqrt(a_weight * a_weight - 1.0);
  const double T = (a_weight < 1.0) ? atan((1.0 - a_weight) / S) / S
                                    : atanh((a_weight - 1.0) / S) / S;
  const GradientType T_grad =
      (a_weight < 1.0)
          ? a_weight_grad * (0.5 + a_weight * T) / (1.0 - a_weight * a_weight)
          : a_weight_grad * (0.5 - a_weight * T) / (a_weight * a_weight - 1.0);
  return std::array<dbl_gradient, 3>(
      {dbl_gradient({(2.0 * w6 - 3.0 * w4 + 31.0 * w2 -
                      (42.0 * w3 + 18.0 * a_weight) * T) *
                         L3 / 12.0,
                     (2.0 * w6_grad - 3.0 * w4_grad + 31.0 * w2_grad -
                      (42.0 * w3_grad + 18.0 * a_weight_grad) * T -
                      (42.0 * w3 + 18.0 * a_weight) * T_grad) *
                             L3 / 12.0 +
                         (2.0 * w6 - 3.0 * w4 + 31.0 * w2 -
                          (42.0 * w3 + 18.0 * a_weight) * T) *
                             L3_grad / 12.0}),
       dbl_gradient(
           {(2.0 * w6 - 9.0 * w4 - 8.0 * w2 + 30.0 * w3 * T) * L3 / 3.0,
            (2.0 * w6_grad - 9.0 * w4_grad - 8.0 * w2_grad +
             30.0 * w3_grad * T + 30.0 * w3 * T_grad) *
                    L3 / 3.0 +
                (2.0 * w6 - 9.0 * w4 - 8.0 * w2 + 30.0 * w3 * T) * L3_grad /
                    3.0}),
       dbl_gradient(
           {(11.0 * w4 + 4.0 * w2 - (12.0 * w5 + 18.0 * w3) * T) * L3 / 6.0,
            (11.0 * w4_grad + 4.0 * w2_grad -
             (12.0 * w5_grad + 18.0 * w3_grad) * T -
             (12.0 * w5 + 18.0 * w3) * T_grad) *
                    L3 / 6.0 +
                (11.0 * w4 + 4.0 * w2 - (12.0 * w5 + 18.0 * w3) * T) * L3_grad /
                    6.0})});
}

inline std::array<double, 12> coeffsV3andC3Exact(const double a_weight) {
  const double w2 = a_weight * a_weight;
  const double w3 = w2 * a_weight;
  const double w4 = w2 * w2;
  const double w5 = w2 * w3;
  const double w6 = w3 * w3;
  const double w7 = w4 * w3;
  const double w8 = w4 * w4;
  const double w9 = w4 * w5;
  const double w10 = w5 * w5;
  const double L = 1.0 / ((a_weight - 1.0) * (a_weight + 1.0));
  const double L3 = L * L * L;
  const double L4 = L3 * L;
  const double L5 = L4 * L;
  const double S = (a_weight < 1.0) ? sqrt(1.0 - a_weight * a_weight)
                                    : sqrt(a_weight * a_weight - 1.0);
  const double T = (a_weight < 1.0) ? atan((1.0 - a_weight) / S) / S
                                    : atanh((a_weight - 1.0) / S) / S;
  return std::array<double, 12>(
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
computeType3Contribution<Volume>(const AlignedParaboloid& a_paraboloid,
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

template <class ReturnType>
inline enable_if_t<has_embedded_gradient<ReturnType>::value, ReturnType>
computeType3Contribution(
    const AlignedParaboloid& a_paraboloid,
    const RationalBezierArcWithGradient<
        PtWithGradient<typename ReturnType::gradient_type>>& a_arc) {
  using gradient_type = typename ReturnType::gradient_type;
  using dbl_gradient = typename std::pair<double, gradient_type>;
  auto volume_with_gradient = ReturnType::fromScalarConstant(0.0);
  const auto& pt_0 = a_arc.start_point().getPt();
  const auto& pt_0_grad = a_arc.start_point().getData();
  const auto& cp = a_arc.control_point().getPt();
  const auto& cp_grad = a_arc.control_point().getData();
  const auto& pt_1 = a_arc.end_point().getPt();
  const auto& pt_1_grad = a_arc.end_point().getData();
  const auto& weight = a_arc.weight();
  assert(weight >= 0.0);
  const auto& weight_grad = a_arc.weight_gradient();
  const auto A = a_paraboloid.a(), B = a_paraboloid.b();
  auto A_grad = gradient_type(0.0), B_grad = gradient_type(0.0);
  A_grad.setGradA(1.0);
  B_grad.setGradB(1.0);
  const double area_proj_triangle =
      0.5 * (pt_0[0] * (pt_1[1] - cp[1]) + pt_1[0] * (cp[1] - pt_0[1]) +
             cp[0] * (pt_0[1] - pt_1[1]));
  const auto area_proj_triangle_grad =
      0.5 * (pt_0_grad[0] * (pt_1[1] - cp[1]) +
             pt_0[0] * (pt_1_grad[1] - cp_grad[1]) +
             pt_1_grad[0] * (cp[1] - pt_0[1]) +
             pt_1[0] * (cp_grad[1] - pt_0_grad[1]) +
             cp_grad[0] * (pt_0[1] - pt_1[1]) +
             cp[0] * (pt_0_grad[1] - pt_1_grad[1]));
  const auto mid1_pt_withgrad = 0.5 * (a_arc.start_point() + a_arc.end_point());
  const auto& mid1_pt = mid1_pt_withgrad.getPt();
  const auto& mid1_grad = mid1_pt_withgrad.getData();
  const auto mid2_pt_withgrad =
      0.5 * (mid1_pt_withgrad + a_arc.control_point());
  const auto& mid2_pt = mid2_pt_withgrad.getPt();
  const auto& mid2_grad = mid2_pt_withgrad.getData();
  auto m0_basis = std::array<double, 3>(
      {area_proj_triangle * signedDistance(mid1_pt, a_paraboloid),
       area_proj_triangle * signedDistance(mid2_pt, a_paraboloid),
       area_proj_triangle * signedDistance(cp, a_paraboloid)});
  auto m0_basis_grad = std::array<gradient_type, 3>(
      {area_proj_triangle_grad * signedDistance(mid1_pt, a_paraboloid) +
           area_proj_triangle *
               (2.0 * A * mid1_grad[0] * mid1_pt[0] +
                A_grad * mid1_pt[0] * mid1_pt[0] +
                2.0 * B * mid1_grad[1] * mid1_pt[1] +
                B_grad * mid1_pt[1] * mid1_pt[1] + mid1_grad[2]),
       area_proj_triangle_grad * signedDistance(mid2_pt, a_paraboloid) +
           area_proj_triangle *
               (2.0 * A * mid2_grad[0] * mid2_pt[0] +
                A_grad * mid2_pt[0] * mid2_pt[0] +
                2.0 * B * mid2_grad[1] * mid2_pt[1] +
                B_grad * mid2_pt[1] * mid2_pt[1] + mid2_grad[2]),
       area_proj_triangle_grad * signedDistance(cp, a_paraboloid) +
           area_proj_triangle *
               (2.0 * A * cp_grad[0] * cp[0] + A_grad * cp[0] * cp[0] +
                2.0 * B * cp_grad[1] * cp[1] + B_grad * cp[1] * cp[1] +
                cp_grad[2])});
  std::array<dbl_gradient, 3> coeffs;
  if (weight < 0.35)  // We use the exact expressions
    coeffs = coeffsV3ExactWithGradient<gradient_type>(weight, weight_grad);
  else if (weight < 1.7)  // We use the 40th order Taylor series (w -> 1)
    coeffs = coeffsV3SeriesOneWithGradient<gradient_type>(weight, weight_grad);
  else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
    coeffs = coeffsV3ExactWithGradient<gradient_type>(weight, weight_grad);
  else  // This is within DBL_EPSILON of the actual value
    coeffs = std::array<dbl_gradient, 3>({dbl_gradient({1.0 / 6.0, 0.0}),
                                          dbl_gradient({2.0 / 3.0, 0.0}),
                                          dbl_gradient({0.0, 0.0})});
  for (size_t i = 0; i < 3; ++i) {
    volume_with_gradient.volume() += coeffs[i].first * m0_basis[i];
    volume_with_gradient.gradient() +=
        coeffs[i].first * m0_basis_grad[i] + coeffs[i].second * m0_basis[i];
  }
  return volume_with_gradient;
}

template <>
inline enable_if_t<!has_embedded_gradient<VolumeMoments>::value, VolumeMoments>
computeType3Contribution<VolumeMoments>(const AlignedParaboloid& a_paraboloid,
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
    coeffs = coeffsV3andC3Exact(weight);
  else if (weight < 1.7)  // We use the 40th order Taylor series (w -> 1)
    coeffs = coeffsV3andC3SeriesOne(weight);
  else if (weight < 1.0e9)  // We use the exact expressions
    coeffs = coeffsV3andC3Exact(weight);
  // else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
  //   coeffs = coeffsV3andC3SeriesInfinity(weight);
  else  // This is within DBL_EPSILON of the actual value
    coeffs = std::array<double, 12>({1.0 / 6.0, 2.0 / 3.0, -1.0 / 140.0,
                                     1.0 / 420.0, -1.0 / 280.0, 1.0 / 14.0, 0.0,
                                     1.0 / 630.0, -1.0 / 1890.0, 1.0 / 3780.0,
                                     -1.0 / 252.0, 1.0 / 18.0});
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
  return moments;
}

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

template <class ReturnType>
inline enable_if_t<has_embedded_gradient<ReturnType>::value, ReturnType>
computeFaceOnlyContribution(
    const AlignedParaboloid& a_paraboloid, const Plane& a_face_plane,
    const PtWithGradient<typename ReturnType::gradient_type>& a_pt_ref) {
  using gradient_type = typename ReturnType::gradient_type;
  auto volume_with_gradient = ReturnType::fromScalarConstant(0.0);
  assert(a_paraboloid.a() * a_paraboloid.b() > 0.0);
  assert(std::fabs(a_face_plane.normal()[2]) > DBL_EPSILON);

  // Compute plane gradient
  const auto& pt_ref = a_pt_ref.getPt();
  const auto& pt_ref_grad = a_pt_ref.getData();
  auto face_normal = a_face_plane.normal();
  auto face_normal_withgrad = PtWithGradient<gradient_type>(
      Pt(face_normal[0], face_normal[1], face_normal[2]));
  auto& face_normal_grad = face_normal_withgrad.getData();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    face_normal_grad[d].setGradRx(0.0);
    face_normal_grad[d].setGradRy(0.0);
    face_normal_grad[d].setGradRz(0.0);
  }
  auto face_distance = face_normal[0] * pt_ref[0] + face_normal[1] * pt_ref[1] +
                       face_normal[2] * pt_ref[2];
  auto face_distance_grad =
      face_normal_grad[0] * pt_ref[0] + face_normal_grad[1] * pt_ref[1] +
      face_normal_grad[2] * pt_ref[2] + face_normal[0] * pt_ref_grad[0] +
      face_normal[1] * pt_ref_grad[1] + face_normal[2] * pt_ref_grad[2];

  const double a = -face_normal[0] / safelyEpsilon(face_normal[2]);
  const double b = -face_normal[1] / safelyEpsilon(face_normal[2]);
  const double c = face_distance / safelyEpsilon(face_normal[2]);
  const auto a_grad =
      -face_normal_grad[0] / safelyEpsilon(face_normal[2]) +
      face_normal[0] / safelyEpsilon(face_normal[2] * face_normal[2]);
  const auto b_grad =
      -face_normal_grad[1] / safelyEpsilon(face_normal[2]) +
      face_normal[1] / safelyEpsilon(face_normal[2] * face_normal[2]);
  const auto c_grad =
      face_distance_grad / safelyEpsilon(face_normal[2]) -
      face_distance / safelyEpsilon(face_normal[2] * face_normal[2]);

  const auto A = a_paraboloid.a(), B = a_paraboloid.b();
  const auto AB = A * B;
  auto A_grad = gradient_type(0.0), B_grad = gradient_type(0.0);
  A_grad.setGradA(1.0);
  B_grad.setGradB(1.0);
  const auto AB_grad = A_grad * B + A * B_grad;
  const double factor = 4.0 * AB * c - A * b * b - B * a * a;
  const auto factor_grad = 4.0 * AB_grad * c + 4.0 * AB * c_grad -
                           A_grad * b * b - 2.0 * A * b_grad * b -
                           B_grad * a * a - 2.0 * B * a_grad * a;
  volume_with_gradient.volume() =
      std::copysign(M_PI * factor * factor / (32.0 * std::pow(AB, 2.5)),
                    -a_face_plane.normal()[2]);
  volume_with_gradient.gradient() =
      std::copysign(1.0, -a_face_plane.normal()[2]) * M_PI *
      (2.0 * factor_grad * factor / (32.0 * std::pow(AB, 2.5)) -
       5.0 * factor * factor * AB_grad / (64.0 * std::pow(AB, 3.5)));
  return volume_with_gradient;
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

  const double x_center = a / (2.0 * a_paraboloid.a());
  const double y_center = b / (2.0 * a_paraboloid.b());
  const double K = 4.0 * a_paraboloid.a() * a_paraboloid.b() * c -
                   a_paraboloid.a() * b * b - a_paraboloid.b() * a * a;

  auto moments = VolumeMoments::fromScalarConstant(0.0);
  moments.volume() = std::copysign(
      M_PI * K * K /
          (32.0 * std::pow(a_paraboloid.a() * a_paraboloid.b(), 2.5)),
      -a_face_plane.normal()[2]);
  moments.centroid()[0] = x_center;
  moments.centroid()[1] = y_center;
  moments.centroid()[2] =
      (8.0 * a_paraboloid.a() * a_paraboloid.b() * c +
       5.0 * a_paraboloid.a() * b * b + 5.0 * a_paraboloid.b() * a * a) /
      (12.0 * a_paraboloid.a() * a_paraboloid.b());
  moments.multiplyByVolume();
  return moments;
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
  auto volume_and_gradient = ReturnType::fromScalarConstant(0.0);
  const auto& pt_0 = a_pt_0.getPt();
  const auto& pt_1 = a_pt_1.getPt();
  const auto& pt_2 = a_pt_2.getPt();
  const auto& pt_0_grad = a_pt_0.getData();
  const auto& pt_1_grad = a_pt_1.getData();
  const auto& pt_2_grad = a_pt_2.getData();
  const double A = a_paraboloid.a(), B = a_paraboloid.b();
  auto A_grad = gradient_type(0.0), B_grad = gradient_type(0.0);
  A_grad.setGradA(1.0);
  B_grad.setGradB(1.0);

  const double triangle_area =
      0.5 * ((pt_1[1] - pt_2[1]) * pt_0[0] + (pt_2[1] - pt_0[1]) * pt_1[0] +
             (pt_0[1] - pt_1[1]) * pt_2[0]);
  const auto triangle_area_grad =
      0.5 *
      ((pt_1_grad[1] - pt_2_grad[1]) * pt_0[0] +
       (pt_2_grad[1] - pt_0_grad[1]) * pt_1[0] +
       (pt_0_grad[1] - pt_1_grad[1]) * pt_2[0] +
       (pt_1[1] - pt_2[1]) * pt_0_grad[0] + (pt_2[1] - pt_0[1]) * pt_1_grad[0] +
       (pt_0[1] - pt_1[1]) * pt_2_grad[0]);
  volume_and_gradient.volume() =
      (-A * (pt_0[0] + pt_1[0]) * (pt_1[0] + pt_2[0]) -
       B * (pt_0[1] + pt_1[1]) * (pt_1[1] + pt_2[1]) - pt_0[2] - 2.0 * pt_1[2] -
       pt_2[2]) *
      triangle_area / 6.0;
  volume_and_gradient.gradient() =
      (-A_grad * (pt_0[0] + pt_1[0]) * (pt_1[0] + pt_2[0]) -
       A * (pt_0_grad[0] + pt_1_grad[0]) * (pt_1[0] + pt_2[0]) -
       A * (pt_0[0] + pt_1[0]) * (pt_1_grad[0] + pt_2_grad[0]) -
       B_grad * (pt_0[1] + pt_1[1]) * (pt_1[1] + pt_2[1]) -
       B * (pt_0_grad[1] + pt_1_grad[1]) * (pt_1[1] + pt_2[1]) -
       B * (pt_0[1] + pt_1[1]) * (pt_1_grad[1] + pt_2_grad[1]) - pt_0_grad[2] -
       2.0 * pt_1_grad[2] - pt_2_grad[2]) *
          triangle_area / 6.0 +
      (-A * (pt_0[0] + pt_1[0]) * (pt_1[0] + pt_2[0]) -
       B * (pt_0[1] + pt_1[1]) * (pt_1[1] + pt_2[1]) - pt_0[2] - 2.0 * pt_1[2] -
       pt_2[2]) *
          triangle_area_grad / 6.0;
  return volume_and_gradient;
}

template <>
inline enable_if_t<!has_embedded_gradient<VolumeMoments>::value, VolumeMoments>
computeTriangleCorrection(const AlignedParaboloid& a_paraboloid,
                          const Pt& a_pt_0, const Pt& a_pt_1,
                          const Pt& a_pt_2) {
  auto moments = VolumeMoments::fromScalarConstant(0.0);
  const double triangle_area = 0.5 * ((a_pt_1[1] - a_pt_2[1]) * a_pt_0[0] +
                                      (a_pt_2[1] - a_pt_0[1]) * a_pt_1[0] +
                                      (a_pt_0[1] - a_pt_1[1]) * a_pt_2[0]);
  moments.volume() =
      (-a_paraboloid.a() * (a_pt_0[0] + a_pt_1[0]) * (a_pt_1[0] + a_pt_2[0]) +
       -a_paraboloid.b() * (a_pt_0[1] + a_pt_1[1]) * (a_pt_1[1] + a_pt_2[1]) -
       a_pt_0[2] - 2.0 * a_pt_1[2] - a_pt_2[2]) *
      triangle_area / 6.0;
  moments.centroid()[0] =
      (6. * a_paraboloid.a() * a_pt_0[0] * a_pt_1[0] * a_pt_2[0] +
       2. * a_paraboloid.b() *
           (a_pt_0[0] * (-2. * (a_pt_1[1] * a_pt_1[1]) + a_pt_1[1] * a_pt_2[1] -
                         2. * (a_pt_2[1] * a_pt_2[1]) +
                         2. * a_pt_0[1] * (a_pt_1[1] + a_pt_2[1])) +
            a_pt_1[0] * (-2. * (a_pt_0[1] * a_pt_0[1]) +
                         2. * (a_pt_1[1] - a_pt_2[1]) * a_pt_2[1] +
                         a_pt_0[1] * (2. * a_pt_1[1] + a_pt_2[1])) +
            a_pt_2[0] * (-2. * (a_pt_0[1] * a_pt_0[1]) +
                         2. * a_pt_1[1] * (-a_pt_1[1] + a_pt_2[1]) +
                         a_pt_0[1] * (a_pt_1[1] + 2. * a_pt_2[1]))) +
       4. * a_pt_0[0] * a_pt_0[2] - a_pt_1[0] * a_pt_0[2] -
       a_pt_2[0] * a_pt_0[2] - a_pt_0[0] * a_pt_1[2] +
       4. * a_pt_1[0] * a_pt_1[2] - a_pt_2[0] * a_pt_1[2] -
       a_pt_0[0] * a_pt_2[2] - a_pt_1[0] * a_pt_2[2] +
       4. * a_pt_2[0] * a_pt_2[2]) *
      triangle_area / 60.;
  moments.centroid()[1] =
      (6. * a_paraboloid.b() * a_pt_0[1] * a_pt_1[1] * a_pt_2[1] +
       2. * a_paraboloid.a() *
           (a_pt_0[1] * (-2. * (a_pt_1[0] * a_pt_1[0]) + a_pt_1[0] * a_pt_2[0] -
                         2. * (a_pt_2[0] * a_pt_2[0]) +
                         2. * a_pt_0[0] * (a_pt_1[0] + a_pt_2[0])) +
            a_pt_1[1] * (-2. * (a_pt_0[0] * a_pt_0[0]) +
                         2. * (a_pt_1[0] - a_pt_2[0]) * a_pt_2[0] +
                         a_pt_0[0] * (2. * a_pt_1[0] + a_pt_2[0])) +
            a_pt_2[1] * (-2. * (a_pt_0[0] * a_pt_0[0]) +
                         2. * a_pt_1[0] * (-a_pt_1[0] + a_pt_2[0]) +
                         a_pt_0[0] * (a_pt_1[0] + 2. * a_pt_2[0]))) +
       4. * a_pt_0[1] * a_pt_0[2] - a_pt_1[1] * a_pt_0[2] -
       a_pt_2[1] * a_pt_0[2] - a_pt_0[1] * a_pt_1[2] +
       4. * a_pt_1[1] * a_pt_1[2] - a_pt_2[1] * a_pt_1[2] -
       a_pt_0[1] * a_pt_2[2] - a_pt_1[1] * a_pt_2[2] +
       4. * a_pt_2[1] * a_pt_2[2]) *
      triangle_area / 60.;
  moments.centroid()[2] =
      ((3. * (a_paraboloid.a() * a_paraboloid.a()) * a_pt_0[0] * a_pt_1[0] *
            a_pt_2[0] * (a_pt_0[0] + a_pt_1[0] + a_pt_2[0]) +
        3. * (a_paraboloid.b() * a_paraboloid.b()) * a_pt_1[1] *
            (a_pt_0[1] * a_pt_0[1] * a_pt_2[1] +
             a_pt_1[1] * a_pt_1[1] * a_pt_2[1] +
             a_pt_0[1] * (a_pt_1[1] * a_pt_1[1] + a_pt_1[1] * a_pt_2[1] +
                          a_pt_2[1] * a_pt_2[1])) -
        3. * a_paraboloid.b() *
            (a_pt_0[1] * a_pt_1[1] * a_pt_0[2] +
             a_pt_1[1] * a_pt_2[1] * a_pt_2[2] +
             a_pt_0[1] * a_pt_2[1] * (a_pt_0[2] + a_pt_2[2])) -
        (9. * (a_pt_0[2] * a_pt_0[2] + a_pt_1[2] * a_pt_1[2] +
               a_pt_1[2] * a_pt_2[2] + a_pt_2[2] * a_pt_2[2] +
               a_pt_0[2] * (a_pt_1[2] + a_pt_2[2]))) /
            2. +
        a_paraboloid.a() *
            (a_paraboloid.b() *
                 (a_pt_2[0] * a_pt_2[0] *
                      (-2. * (a_pt_0[1] * a_pt_0[1]) + a_pt_0[1] * a_pt_1[1] -
                       2. * (a_pt_1[1] * a_pt_1[1])) +
                  a_pt_0[0] * a_pt_2[0] * (2. * a_pt_0[1] + a_pt_1[1]) *
                      (a_pt_1[1] + 2. * a_pt_2[1]) +
                  a_pt_0[0] * a_pt_0[0] *
                      (-2. * (a_pt_1[1] * a_pt_1[1]) + a_pt_1[1] * a_pt_2[1] -
                       2. * (a_pt_2[1] * a_pt_2[1])) +
                  a_pt_1[0] * a_pt_1[0] *
                      (-2. * (a_pt_0[1] * a_pt_0[1]) +
                       (3. * a_pt_1[1] - 2. * a_pt_2[1]) * a_pt_2[1] +
                       a_pt_0[1] * (3. * a_pt_1[1] + a_pt_2[1])) +
                  a_pt_1[0] * (a_pt_0[0] * (2. * a_pt_0[1] + a_pt_2[1]) *
                                   (2. * a_pt_1[1] + a_pt_2[1]) +
                               a_pt_2[0] * (a_pt_0[1] + 2. * a_pt_1[1]) *
                                   (a_pt_0[1] + 2. * a_pt_2[1]))) -
             3. * (a_pt_0[0] * a_pt_1[0] * (a_pt_0[2] + a_pt_1[2]) +
                   a_pt_0[0] * a_pt_2[0] * (a_pt_0[2] + a_pt_2[2]) +
                   a_pt_1[0] * a_pt_2[0] * (a_pt_1[2] + a_pt_2[2]))))) *
      triangle_area / 90.;
  return moments;
}

}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_MOMENT_CONTRIBUTIONS_TPP_
