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
template <class ScalarType>
inline ScalarType signedDistance(
    const PtBase<ScalarType>& a_pt,
    const AlignedParaboloidBase<ScalarType>& a_paraboloid) {
  return a_paraboloid.a() * a_pt[0] * a_pt[0] +
         a_paraboloid.b() * a_pt[1] * a_pt[1] + a_pt[2];
}

/******************************************************************************/
/*********************** First moment contribution ****************************/
/******************************************************************************/
/* This compute the first contribution to the moments (arising from the
 * integration of the face plane primitives on the poligonized clipped faces) */
template <class ReturnType, class ScalarType>
ReturnType computeType1Contribution(const PtBase<ScalarType>& a_ref_pt,
                                    const PtBase<ScalarType>& a_pt_0,
                                    const PtBase<ScalarType>& a_pt_1) {
  if constexpr (std::is_same_v<ReturnType, VolumeBase<double>> ||
                std::is_same_v<ReturnType, VolumeBase<Quad_t>>) {
    return (a_ref_pt[2] + a_pt_0[2] + a_pt_1[2]) /
           static_cast<ScalarType>(6.0) *
           ((a_pt_0[0] - a_ref_pt[0]) * (a_pt_1[1] - a_ref_pt[1]) -
            (a_pt_1[0] - a_ref_pt[0]) * (a_pt_0[1] - a_ref_pt[1]));
  } else if constexpr (std::is_same_v<ReturnType, VolumeMomentsBase<double>> ||
                       std::is_same_v<ReturnType, VolumeMomentsBase<Quad_t>>) {
    /* Defining constants and types */
    const ScalarType ZERO = static_cast<ScalarType>(0);
    const ScalarType ONE = static_cast<ScalarType>(1);
    const ScalarType TWO = static_cast<ScalarType>(2);
    const ScalarType THREE = static_cast<ScalarType>(3);
    const ScalarType ONETWELVTH = ONE / static_cast<ScalarType>(12);

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ZERO);
    const ScalarType triangle_area =
        ((a_pt_0[0] - a_ref_pt[0]) * (a_pt_1[1] - a_ref_pt[1]) -
         (a_pt_1[0] - a_ref_pt[0]) * (a_pt_0[1] - a_ref_pt[1])) /
        TWO;
    moments.volume() =
        triangle_area * (a_ref_pt[2] + a_pt_0[2] + a_pt_1[2]) / THREE;
    moments.centroid()[0] =
        triangle_area *
        (a_pt_0[0] * (TWO * a_pt_0[2] + a_pt_1[2] + a_ref_pt[2]) +
         a_pt_1[0] * (a_pt_0[2] + TWO * a_pt_1[2] + a_ref_pt[2]) +
         a_ref_pt[0] * (a_pt_0[2] + a_pt_1[2] + TWO * a_ref_pt[2])) *
        ONETWELVTH;
    moments.centroid()[1] =
        triangle_area *
        (a_pt_0[1] * (TWO * a_pt_0[2] + a_pt_1[2] + a_ref_pt[2]) +
         a_pt_1[1] * (a_pt_0[2] + TWO * a_pt_1[2] + a_ref_pt[2]) +
         a_ref_pt[1] * (a_pt_0[2] + a_pt_1[2] + TWO * a_ref_pt[2])) *
        ONETWELVTH;
    moments.centroid()[2] =
        triangle_area *
        (a_pt_0[2] * a_pt_0[2] + a_pt_1[2] * a_pt_1[2] +
         a_ref_pt[2] * a_ref_pt[2] + a_pt_1[2] * a_ref_pt[2] +
         a_pt_0[2] * a_pt_1[2] + a_pt_0[2] * a_ref_pt[2]) *
        ONETWELVTH;
    return moments;
  }
}

/* This is the version of computeType1Contribution for returning the zeroth
 * moment only or the zeroth and first moments with their gradient */
template <class ReturnType, class ScalarType>
inline ReturnType computeType1ContributionWithGradient(
    const PtWithGradientBase<typename ReturnType::gradient_type, ScalarType>&
        a_ref_pt,
    const PtWithGradientBase<typename ReturnType::gradient_type, ScalarType>&
        a_pt_0,
    const PtWithGradientBase<typename ReturnType::gradient_type, ScalarType>&
        a_pt_1) {
  /* Defining constants and types */
  using gradient_type = typename ReturnType::gradient_type;
  using scalar_withgrad_type = ScalarWithGradient<gradient_type>;
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType THREE = static_cast<ScalarType>(3);
  const ScalarType ONETWELVTH = ONE / static_cast<ScalarType>(12);

  /* Function */
  auto moments_with_gradient = ReturnType::fromScalarConstant(ZERO);
  std::array<scalar_withgrad_type, 3> pt_0, pt_1, ref_pt;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_0[d] = scalar_withgrad_type(a_pt_0.getPt()[d], a_pt_0.getData()[d]);
    pt_1[d] = scalar_withgrad_type(a_pt_1.getPt()[d], a_pt_1.getData()[d]);
    ref_pt[d] =
        scalar_withgrad_type(a_ref_pt.getPt()[d], a_ref_pt.getData()[d]);
  }
  const auto triangle_area = ((pt_0[0] - ref_pt[0]) * (pt_1[1] - ref_pt[1]) -
                              (pt_1[0] - ref_pt[0]) * (pt_0[1] - ref_pt[1])) /
                             TWO;
  const auto volume = triangle_area * (ref_pt[2] + pt_0[2] + pt_1[2]) / THREE;
  moments_with_gradient.volume() = volume.value();
  moments_with_gradient.volume_gradient() = volume.gradient();
  if constexpr (std::is_same<typename ReturnType::moment_type,
                             VolumeMoments>::value) {
    PtBase<ScalarType>& centroid = moments_with_gradient.centroid().getPt();
    auto& centroid_gradient = moments_with_gradient.centroid().getData();
    const auto m1x = triangle_area *
                     (pt_0[0] * (TWO * pt_0[2] + pt_1[2] + ref_pt[2]) +
                      pt_1[0] * (pt_0[2] + TWO * pt_1[2] + ref_pt[2]) +
                      ref_pt[0] * (pt_0[2] + pt_1[2] + TWO * ref_pt[2])) *
                     ONETWELVTH;
    const auto m1y = triangle_area *
                     (pt_0[1] * (TWO * pt_0[2] + pt_1[2] + ref_pt[2]) +
                      pt_1[1] * (pt_0[2] + TWO * pt_1[2] + ref_pt[2]) +
                      ref_pt[1] * (pt_0[2] + pt_1[2] + TWO * ref_pt[2])) *
                     ONETWELVTH;
    const auto m1z =
        triangle_area *
        (pt_0[2] * pt_0[2] + pt_1[2] * pt_1[2] + ref_pt[2] * ref_pt[2] +
         pt_1[2] * ref_pt[2] + pt_0[2] * pt_1[2] + pt_0[2] * ref_pt[2]) *
        ONETWELVTH;
    centroid[0] = m1x.value();
    centroid[1] = m1y.value();
    centroid[2] = m1z.value();
    centroid_gradient[0] = m1x.gradient();
    centroid_gradient[1] = m1y.gradient();
    centroid_gradient[2] = m1z.gradient();
  }
  return moments_with_gradient;
}

template <class ReturnType, class ScalarType>
ReturnType computeType2Contribution(
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid,
    const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1) {
  if constexpr (std::is_same_v<ReturnType, VolumeBase<double>> ||
                std::is_same_v<ReturnType, VolumeBase<Quad_t>>) {
    const ScalarType ONETWELVTH =
        static_cast<ScalarType>(1) / static_cast<ScalarType>(12);
    return ONETWELVTH * (a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) *
           (-a_pt_0[2] - a_pt_1[2] +
            a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] +
            a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1]);
  } else if constexpr (std::is_same_v<ReturnType, VolumeMomentsBase<double>> ||
                       std::is_same_v<ReturnType, VolumeMomentsBase<Quad_t>>) {
    /* Defining constants and types */
    const ScalarType ZERO = static_cast<ScalarType>(0);
    const ScalarType ONE = static_cast<ScalarType>(1);
    const ScalarType TWO = static_cast<ScalarType>(2);
    const ScalarType THREE = static_cast<ScalarType>(3);
    const ScalarType ONETWELVTH = ONE / static_cast<ScalarType>(12);
    const ScalarType ONE60TH = ONE / static_cast<ScalarType>(60);
    const ScalarType ONE180TH = ONE / static_cast<ScalarType>(180);

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ZERO);
    moments.volume() = (a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) *
                       ONETWELVTH *
                       (-a_pt_0[2] - a_pt_1[2] +
                        a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] +
                        a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1]);
    moments.centroid()[0] =
        (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) *
        (TWO * a_aligned_paraboloid.b() * (a_pt_0[1] - a_pt_1[1]) *
             (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) +
         THREE * (a_pt_0[0] + a_pt_1[0]) * (a_pt_0[2] + a_pt_1[2])) *
        ONE60TH;
    moments.centroid()[1] =
        (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) *
        (TWO * a_aligned_paraboloid.a() * (a_pt_0[0] - a_pt_1[0]) *
             (a_pt_1[1] * a_pt_0[0] - a_pt_0[1] * a_pt_1[0]) +
         THREE * (a_pt_0[1] + a_pt_1[1]) * (a_pt_0[2] + a_pt_1[2])) *
        ONE60TH;
    moments.centroid()[2] =
        ((a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) *
         (TWO * a_aligned_paraboloid.a() * a_aligned_paraboloid.b() *
              ((a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) *
               (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1])) +
          THREE * a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] *
              (a_pt_0[2] + a_pt_1[2]) +
          THREE * a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1] *
              (a_pt_0[2] + a_pt_1[2]) -
          THREE * (a_pt_0[2] * a_pt_0[2] + a_pt_0[2] * a_pt_1[2] +
                   a_pt_1[2] * a_pt_1[2]))) *
        ONE180TH;
    return moments;
  }
}

template <class ReturnType, class ScalarType>
inline ReturnType computeType2ContributionWithGradient(
    const AlignedParaboloid& a_aligned_paraboloid,
    const PtWithGradientBase<typename ReturnType::gradient_type, ScalarType>&
        a_pt_0,
    const PtWithGradientBase<typename ReturnType::gradient_type, ScalarType>&
        a_pt_1) {
  /* Defining constants and types */
  using gradient_type = typename ReturnType::gradient_type;
  using scalar_container_type = ScalarWithGradient<gradient_type>;
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType THREE = static_cast<ScalarType>(3);
  const ScalarType ONETWELVTH = ONE / static_cast<ScalarType>(12);
  const ScalarType ONE60TH = ONE / static_cast<ScalarType>(60);
  const ScalarType ONE180TH = ONE / static_cast<ScalarType>(180);

  /* Function */
  auto moments_with_gradient = ReturnType::fromScalarConstant(ZERO);
  std::array<scalar_container_type, 3> pt_0, pt_1;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_0[d] = scalar_container_type(a_pt_0.getPt()[d], a_pt_0.getData()[d]);
    pt_1[d] = scalar_container_type(a_pt_1.getPt()[d], a_pt_1.getData()[d]);
  }
  auto A_grad = gradient_type(ZERO), B_grad = gradient_type(ZERO);
  A_grad.setGradA(ONE);
  B_grad.setGradB(ONE);
  const auto A = scalar_container_type(a_aligned_paraboloid.a(), A_grad),
             B = scalar_container_type(a_aligned_paraboloid.b(), B_grad);
  const auto volume =
      (pt_0[0] * pt_1[1] - pt_1[0] * pt_0[1]) * ONETWELVTH *
      (-pt_0[2] - pt_1[2] + A * pt_0[0] * pt_1[0] + B * pt_0[1] * pt_1[1]);
  moments_with_gradient.volume() = volume.value();
  moments_with_gradient.volume_gradient() = volume.gradient();
  if constexpr (std::is_same<typename ReturnType::moment_type,
                             VolumeMoments>::value) {
    PtBase<ScalarType>& centroid = moments_with_gradient.centroid().getPt();
    auto& centroid_gradient = moments_with_gradient.centroid().getData();
    const auto m1x = (pt_1[0] * pt_0[1] - pt_0[0] * pt_1[1]) *
                     (TWO * B * (pt_0[1] - pt_1[1]) *
                          (pt_1[0] * pt_0[1] - pt_0[0] * pt_1[1]) +
                      THREE * (pt_0[0] + pt_1[0]) * (pt_0[2] + pt_1[2])) *
                     ONE60TH;
    const auto m1y = (pt_1[0] * pt_0[1] - pt_0[0] * pt_1[1]) *
                     (TWO * A * (pt_0[0] - pt_1[0]) *
                          (pt_1[1] * pt_0[0] - pt_0[1] * pt_1[0]) +
                      THREE * (pt_0[1] + pt_1[1]) * (pt_0[2] + pt_1[2])) *
                     ONE60TH;
    const auto m1z = ((pt_0[0] * pt_1[1] - pt_1[0] * pt_0[1]) *
                      (TWO * A * B *
                           ((pt_1[0] * pt_0[1] - pt_0[0] * pt_1[1]) *
                            (pt_1[0] * pt_0[1] - pt_0[0] * pt_1[1])) +
                       THREE * A * pt_0[0] * pt_1[0] * (pt_0[2] + pt_1[2]) +
                       THREE * B * pt_0[1] * pt_1[1] * (pt_0[2] + pt_1[2]) -
                       THREE * (pt_0[2] * pt_0[2] + pt_0[2] * pt_1[2] +
                                pt_1[2] * pt_1[2]))) *
                     ONE180TH;
    centroid[0] = m1x.value();
    centroid[1] = m1y.value();
    centroid[2] = m1z.value();
    centroid_gradient[0] = m1x.gradient();
    centroid_gradient[1] = m1y.gradient();
    centroid_gradient[2] = m1z.gradient();
  }
  return moments_with_gradient;
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 3> coeffsV3SeriesOne(
    const ContainerType& a_weight) {
  std::array<ContainerType, 3> coeffs;
  coeffs.fill(ContainerType(static_cast<ScalarType>(0)));
  ContainerType x(static_cast<ScalarType>(1));
  UnsignedIndex_t i = 0;
  ScalarType max_diff;
  while (i <= 40) {
    max_diff = static_cast<ScalarType>(0);
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      ContainerType add_to_coeff = static_cast<ScalarType>(v3Series[i][j]) * x;
      coeffs[j] += add_to_coeff;
      if constexpr (has_embedded_gradient<ContainerType>::value) {
        max_diff = maximum(max_diff, fabs(add_to_coeff.value()));
      } else {
        max_diff = maximum(max_diff, fabs(add_to_coeff));
      }
    }
    if (max_diff < static_cast<ScalarType>(10.0 * DBL_EPSILON)) {
      break;
    }
    x *= a_weight - ContainerType(static_cast<ScalarType>(1));
    i++;
  }
  return coeffs;
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 12> coeffsV3andC3SeriesOne(
    const ContainerType& a_weight) {
  std::array<ContainerType, 12> coeffs;
  coeffs.fill(ContainerType(static_cast<ScalarType>(0)));
  ContainerType x(static_cast<ScalarType>(1));
  UnsignedIndex_t i = 0;
  ScalarType max_diff;
  while (i <= 40) {
    max_diff = static_cast<ScalarType>(0);
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      ContainerType add_to_coeff = static_cast<ScalarType>(v3Series[i][j]) * x;
      coeffs[j] += add_to_coeff;
      if constexpr (has_embedded_gradient<ContainerType>::value) {
        max_diff = maximum(max_diff, fabs(add_to_coeff.value()));
      } else {
        max_diff = maximum(max_diff, fabs(add_to_coeff));
      }
    }
    for (UnsignedIndex_t j = 0; j < 4; ++j) {
      ContainerType add_to_coeff = static_cast<ScalarType>(cx3Series[i][j]) * x;
      coeffs[3 + j] += add_to_coeff;
      if constexpr (has_embedded_gradient<ContainerType>::value) {
        max_diff = maximum(max_diff, fabs(add_to_coeff.value()));
      } else {
        max_diff = maximum(max_diff, fabs(add_to_coeff));
      }
    }
    for (UnsignedIndex_t j = 0; j < 5; ++j) {
      ContainerType add_to_coeff = static_cast<ScalarType>(cz3Series[i][j]) * x;
      coeffs[7 + j] += add_to_coeff;
      if constexpr (has_embedded_gradient<ContainerType>::value) {
        max_diff = maximum(max_diff, fabs(add_to_coeff.value()));
      } else {
        max_diff = maximum(max_diff, fabs(add_to_coeff));
      }
    }
    if (max_diff < static_cast<ScalarType>(10.0 * DBL_EPSILON)) {
      break;
    }
    x *= a_weight - ContainerType(static_cast<ScalarType>(1));
    i++;
  }

  return coeffs;
}

// template <class ScalarType>
// inline std::array<ScalarType, 3> coeffsV3SeriesInfinity(
//     const ScalarType& a_weight) {
//   const auto wm2 = 1.0 / (a_weight * a_weight);
//   const auto wm4 = wm2 * wm2;
//   const auto wm6 = wm4 * wm2;
//   const auto wm8 = wm4 * wm4;
//   const auto ln2plnw =
//       ScalarType(std::log(2.0)) + LogMoments<ScalarType>(a_weight);
//   return std::array<ScalarType, 3>(
//       {1.0 / 6.0 + wm2 / 4.0 + wm4 * (17.0 / 6.0 - ln2plnw * 7.0 / 4.0) +
//            wm6 * (401.0 / 48.0 - ln2plnw * 55.0 / 8.0) +
//            wm8 * (2225.0 / 128.0 - ln2plnw * 525.0 / 32.0),
//        2.0 / 3.0 - wm2 + wm4 * (-23.0 / 3.0 + 5.0 * ln2plnw) +
//            wm6 * (-247.0 / 12.0 + 35.0 * ln2plnw / 2.0) +
//            wm8 * (-1307.0 / 32.0 + 315.0 * ln2plnw / 8.0),
//        wm2 * (11.0 / 6.0 - ln2plnw) + wm4 * (77.0 / 12.0 - 5.0 * ln2plnw) +
//            wm6 * (459.0 / 32.0 - 105.0 * ln2plnw / 8.0) +
//            wm8 * (2509.0 / 96.0 - 105.0 * ln2plnw / 4.0)});
// }

// template <class ScalarType>
// inline std::array<ScalarType, 12> coeffsV3andC3SeriesInfinity(
//     const double a_weight) {
//   std::cout << "Not implemented yet" << std::endl;
//   std::exit(-1);
//   return std::array<ScalarType, 12>();
// }

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 3> coeffsV3Exact(
    const ContainerType& a_weight) {
  /* Defining constants and types */
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType THREE = static_cast<ScalarType>(3);
  const ScalarType FOUR = static_cast<ScalarType>(4);
  const ScalarType SIX = static_cast<ScalarType>(6);

  /* Function */
  const auto w2 = a_weight * a_weight;
  const auto w3 = w2 * a_weight;
  const auto w4 = w2 * w2;
  const auto w5 = w2 * w3;
  const auto w6 = w3 * w3;
  const auto L = ContainerType(ONE) / ((a_weight - ContainerType(ONE)) *
                                       (a_weight + ContainerType(ONE)));
  const auto L3 = L * L * L;
  const auto S =
      (a_weight < ContainerType(ONE))
          ? SqrtMoments<ContainerType>(ContainerType(ONE) - a_weight * a_weight)
          : SqrtMoments<ContainerType>(a_weight * a_weight -
                                       ContainerType(ONE));
  const auto T =
      (a_weight < ContainerType(ONE))
          ? ArctanMoments<ContainerType>((ContainerType(ONE) - a_weight) / S) /
                S
          : ArctanhMoments<ContainerType>((a_weight - ContainerType(ONE)) / S) /
                S;
  return std::array<ContainerType, 3>(
      {(TWO * w6 - THREE * w4 + static_cast<ScalarType>(31) * w2 -
        (static_cast<ScalarType>(42) * w3 +
         static_cast<ScalarType>(18) * a_weight) *
            T) *
           L3 / static_cast<ScalarType>(12),
       (TWO * w6 - static_cast<ScalarType>(9) * w4 -
        static_cast<ScalarType>(8) * w2 +
        static_cast<ScalarType>(30) * w3 * T) *
           L3 / THREE,
       (static_cast<ScalarType>(11) * w4 + FOUR * w2 -
        (static_cast<ScalarType>(12) * w5 + static_cast<ScalarType>(18) * w3) *
            T) *
           L3 / SIX});
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 12> coeffsV3andC3Exact(
    const ContainerType& a_weight) {
  /* Defining constants and types */
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType THREE = static_cast<ScalarType>(3);
  const ScalarType FOUR = static_cast<ScalarType>(4);
  const ScalarType FIVE = static_cast<ScalarType>(5);
  const ScalarType SIX = static_cast<ScalarType>(6);

  /* Function */
  const auto w2 = a_weight * a_weight;
  const auto w3 = w2 * a_weight;
  const auto w4 = w2 * w2;
  const auto w5 = w2 * w3;
  const auto w6 = w3 * w3;
  const auto w7 = w4 * w3;
  const auto w8 = w4 * w4;
  const auto w9 = w4 * w5;
  const auto w10 = w5 * w5;
  const auto L = ContainerType(ONE) / ((a_weight - ContainerType(ONE)) *
                                       (a_weight + ContainerType(ONE)));
  const auto L3 = L * L * L;
  const auto L4 = L3 * L;
  const auto L5 = L4 * L;
  const auto S =
      (a_weight < ContainerType(ONE))
          ? SqrtMoments<ContainerType>(ContainerType(ONE) - a_weight * a_weight)
          : SqrtMoments<ContainerType>(a_weight * a_weight -
                                       ContainerType(ONE));
  const auto T =
      (a_weight < ContainerType(ONE))
          ? ArctanMoments<ContainerType>((ContainerType(ONE) - a_weight) / S) /
                S
          : ArctanhMoments<ContainerType>((a_weight - ContainerType(ONE)) / S) /
                S;
  return std::array<ContainerType, 12>(
      {L3 *
           (TWO * w6 - THREE * w4 + static_cast<ScalarType>(31) * w2 -
            (static_cast<ScalarType>(42) * w3 +
             static_cast<ScalarType>(18) * a_weight) *
                T) /
           static_cast<ScalarType>(12),
       L3 *
           (TWO * w6 - static_cast<ScalarType>(9) * w4 -
            static_cast<ScalarType>(8) * w2 +
            static_cast<ScalarType>(30) * w3 * T) /
           THREE,
       L3 *
           (static_cast<ScalarType>(11) * w4 + FOUR * w2 -
            (static_cast<ScalarType>(12) * w5 +
             static_cast<ScalarType>(18) * w3) *
                T) /
           SIX,
       L4 * ((-T * a_weight) / static_cast<ScalarType>(32) +
             (static_cast<ScalarType>(93) * (w2)) /
                 static_cast<ScalarType>(2240) -
             (static_cast<ScalarType>(163) * (w4)) /
                 static_cast<ScalarType>(3360) +
             (FIVE * (w6)) / static_cast<ScalarType>(168) -
             (w8) / static_cast<ScalarType>(140)),
       L4 * ((w2) / static_cast<ScalarType>(70) +
             (-T * (w3)) / static_cast<ScalarType>(16) +
             (static_cast<ScalarType>(29) * (w4)) /
                 static_cast<ScalarType>(1120) -
             (static_cast<ScalarType>(19) * (w6)) /
                 static_cast<ScalarType>(1680) +
             (w8) / static_cast<ScalarType>(420)),
       -L4 * ((w2) / static_cast<ScalarType>(210) -
              (w4) / static_cast<ScalarType>(21) -
              (-T * (w5)) / static_cast<ScalarType>(8) -
              (static_cast<ScalarType>(13) * (w6)) /
                  static_cast<ScalarType>(560) +
              (w8) / static_cast<ScalarType>(280)),
       L4 * ((w2) / static_cast<ScalarType>(35) -
             (static_cast<ScalarType>(16) * (w4)) /
                 static_cast<ScalarType>(105) +
             (static_cast<ScalarType>(58) * (w6)) /
                 static_cast<ScalarType>(105) -
             T * (w7) + (w8) / static_cast<ScalarType>(14)),
       L5 * ((-T * a_weight) / static_cast<ScalarType>(128) +
             (static_cast<ScalarType>(193) * (w2)) /
                 static_cast<ScalarType>(16128) -
             (static_cast<ScalarType>(149) * (w4)) /
                 static_cast<ScalarType>(8064) +
             (static_cast<ScalarType>(19) * (w6)) /
                 static_cast<ScalarType>(1120) -
             (static_cast<ScalarType>(41) * (w8)) /
                 static_cast<ScalarType>(5040) +
             (w10) / static_cast<ScalarType>(630)),
       L5 * ((FOUR * (w2)) / static_cast<ScalarType>(945) +
             (-T * (w3)) / static_cast<ScalarType>(48) +
             (static_cast<ScalarType>(65) * (w4)) /
                 static_cast<ScalarType>(6048) -
             (w6) / static_cast<ScalarType>(144) +
             (static_cast<ScalarType>(11) * (w8)) /
                 static_cast<ScalarType>(3780) -
             (w10) / static_cast<ScalarType>(1890)),
       -L5 * ((w2) / static_cast<ScalarType>(1890) -
              (static_cast<ScalarType>(13) * (w4)) /
                  static_cast<ScalarType>(1890) -
              (-T * (w5)) / static_cast<ScalarType>(48) -
              (static_cast<ScalarType>(11) * (w6)) /
                  static_cast<ScalarType>(2016) +
              (static_cast<ScalarType>(5) * (w8)) /
                  static_cast<ScalarType>(3024) -
              (w10) / static_cast<ScalarType>(3780)),
       L5 * ((w2) / static_cast<ScalarType>(315) -
             (w4) / static_cast<ScalarType>(45) +
             (static_cast<ScalarType>(4) * (w6)) / static_cast<ScalarType>(35) +
             (-T * (w7)) / static_cast<ScalarType>(4) +
             (static_cast<ScalarType>(17) * (w8)) /
                 static_cast<ScalarType>(504) -
             (w10) / static_cast<ScalarType>(252)),
       -L5 * ((w2) / static_cast<ScalarType>(63) -
              (static_cast<ScalarType>(29) * (w4)) /
                  static_cast<ScalarType>(315) +
              (static_cast<ScalarType>(26) * (w6)) /
                  static_cast<ScalarType>(105) -
              (static_cast<ScalarType>(194) * (w8)) /
                  static_cast<ScalarType>(315) +
              T * (w9) - (w10) / static_cast<ScalarType>(18))});
}

template <class ReturnType, class ScalarType>
ReturnType computeType3Contribution(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const RationalBezierArcBase<ScalarType>& a_arc) {
  if constexpr (std::is_same_v<ReturnType, VolumeBase<double>> ||
                std::is_same_v<ReturnType, VolumeBase<Quad_t>>) {
    /* Defining constants and types */
    const ScalarType ZERO = static_cast<ScalarType>(0);
    const ScalarType ONE = static_cast<ScalarType>(1);
    const ScalarType TWO = static_cast<ScalarType>(2);
    const ScalarType THREE = static_cast<ScalarType>(3);
    const ScalarType FOUR = static_cast<ScalarType>(4);
    const ScalarType FIVE = static_cast<ScalarType>(5);
    const ScalarType SIX = static_cast<ScalarType>(6);
    const ScalarType HALF = ONE / TWO;
    const ScalarType QUARTER = ONE / FOUR;

    const auto& pt_0 = a_arc.start_point();
    const auto& cp = a_arc.control_point();
    const auto& pt_1 = a_arc.end_point();
    const auto& weight = a_arc.weight();
    const ScalarType area_proj_triangle =
        HALF * (pt_0[0] * (pt_1[1] - cp[1]) + pt_1[0] * (cp[1] - pt_0[1]) +
                cp[0] * (pt_0[1] - pt_1[1]));
    assert(weight >= ZERO);
    std::array<ScalarType, 3> coeffs;
    if (weight < static_cast<ScalarType>(0.35))  // We use the exact expressions
      coeffs = coeffsV3Exact<ScalarType, ScalarType>(weight);
    else if (weight < static_cast<ScalarType>(
                          1.7))  // We use the 40th order Taylor series (w -> 1)
      coeffs = coeffsV3SeriesOne<ScalarType, ScalarType>(weight);
    else if (weight < static_cast<ScalarType>(
                          1.0e9))  // We use the series expansion (w -> infty)
      coeffs = coeffsV3Exact<ScalarType, ScalarType>(weight);
    else  // This is within EPSILON of the actual value
      coeffs = std::array<ScalarType, 3>({ONE / SIX, TWO / THREE, ZERO});

    return area_proj_triangle *
           (coeffs[0] *
                signedDistance<ScalarType>(HALF * (pt_0 + pt_1), a_paraboloid) +
            coeffs[1] * signedDistance<ScalarType>(
                            QUARTER * (pt_0 + pt_1) + HALF * cp, a_paraboloid) +
            coeffs[2] * signedDistance<ScalarType>(cp, a_paraboloid));
  } else if constexpr (std::is_same_v<ReturnType, VolumeMomentsBase<double>> ||
                       std::is_same_v<ReturnType, VolumeMomentsBase<Quad_t>>) {
    /* Defining constants and types */
    const ScalarType ZERO = static_cast<ScalarType>(0);
    const ScalarType ONE = static_cast<ScalarType>(1);
    const ScalarType TWO = static_cast<ScalarType>(2);
    const ScalarType THREE = static_cast<ScalarType>(3);
    const ScalarType FOUR = static_cast<ScalarType>(4);
    const ScalarType FIVE = static_cast<ScalarType>(5);
    const ScalarType SIX = static_cast<ScalarType>(6);
    const ScalarType HALF = ONE / TWO;
    const ScalarType QUARTER = ONE / FOUR;

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ZERO);
    const auto& pt_0 = a_arc.start_point();
    const auto& cp = a_arc.control_point();
    const auto& pt_1 = a_arc.end_point();
    const auto weight = a_arc.weight();
    const auto A = a_paraboloid.a(), B = a_paraboloid.b();
    const auto X0 = pt_0[0], X1 = cp[0], X2 = pt_1[0];
    const auto Y0 = pt_0[1], Y1 = cp[1], Y2 = pt_1[1];
    const auto Z0 = pt_0[2], Z1p = cp[2], Z2 = pt_1[2];
    const ScalarType Z1P = -A * X1 * X1 - B * Y1 * Y1;
    const ScalarType AA = A * A, BB = B * B, AB = A * B;
    const ScalarType X00 = X0 * X0, X11 = X1 * X1, X22 = X2 * X2;
    const ScalarType X000 = X00 * X0, X111 = X11 * X1, X222 = X22 * X2;
    const ScalarType X0000 = X00 * X00, X1111 = X11 * X11, X2222 = X22 * X22;
    const ScalarType Y00 = Y0 * Y0, Y11 = Y1 * Y1, Y22 = Y2 * Y2;
    const ScalarType Y000 = Y00 * Y0, Y111 = Y11 * Y1, Y222 = Y22 * Y2;
    const ScalarType Y0000 = Y00 * Y00, Y1111 = Y11 * Y11, Y2222 = Y22 * Y22;
    const ScalarType Z00 = Z0 * Z0, Z22 = Z2 * Z2;
    const ScalarType Z1p1p = Z1p * Z1p, Z1P1P = Z1P * Z1P;
    const ScalarType X02 = X0 * X2, X12 = X1 * X2, X01 = X0 * X1;
    const ScalarType Y02 = Y0 * Y2, Y12 = Y1 * Y2, Y01 = Y0 * Y1;
    const ScalarType Z02 = Z0 * Z2, Z01p = Z0 * Z1p, Z01P = Z0 * Z1P;
    const ScalarType Z1p1P = Z1p * Z1P, Z1p2 = Z1p * Z2, Z1P2 = Z1P * Z2;
    const ScalarType X0Z0 = X0 * Z0, X0Z1p = X0 * Z1p, X0Z1P = X0 * Z1P,
                     X0Z2 = X0 * Z2;
    const ScalarType X1Z0 = X1 * Z0, X1Z1p = X1 * Z1p, X1Z1P = X1 * Z1P,
                     X1Z2 = X1 * Z2;
    const ScalarType X2Z0 = X2 * Z0, X2Z1p = X2 * Z1p, X2Z1P = X2 * Z1P,
                     X2Z2 = X2 * Z2;
    const ScalarType Y0Z0 = Y0 * Z0, Y0Z1p = Y0 * Z1p, Y0Z1P = Y0 * Z1P,
                     Y0Z2 = Y0 * Z2;
    const ScalarType Y1Z0 = Y1 * Z0, Y1Z1p = Y1 * Z1p, Y1Z1P = Y1 * Z1P,
                     Y1Z2 = Y1 * Z2;
    const ScalarType Y2Z0 = Y2 * Z0, Y2Z1p = Y2 * Z1p, Y2Z1P = Y2 * Z1P,
                     Y2Z2 = Y2 * Z2;
    const ScalarType area_proj_triangle =
        HALF * (X0 * (Y2 - Y1) + X2 * (Y1 - Y0) + X1 * (Y0 - Y2));
    assert(weight >= ZERO);
    // Compute coefficients (functions of the weight)
    std::array<ScalarType, 12> coeffs;
    if (weight < static_cast<ScalarType>(0.35))  // We use the exact expressions
    {
      coeffs = coeffsV3andC3Exact<ScalarType, ScalarType>(weight);
    } else if (weight <
               static_cast<ScalarType>(
                   1.7))  // We use the 40th order Taylor series (w -> 1)
    {
      coeffs = coeffsV3andC3SeriesOne<ScalarType, ScalarType>(weight);
    } else if (weight <
               static_cast<ScalarType>(1.0e9))  // We use the exact expressions
    {
      coeffs = coeffsV3andC3Exact<ScalarType, ScalarType>(weight);
    }
    // else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
    //   coeffs = coeffsV3andC3SeriesInfinity(weight);
    else  // This is within EPSILON of the actual value
    {
      coeffs = std::array<ScalarType, 12>(
          {ScalarType(ONE / SIX), ScalarType(TWO / THREE),
           ScalarType(-ONE / static_cast<ScalarType>(140)),
           ScalarType(ONE / static_cast<ScalarType>(420)),
           ScalarType(-ONE / static_cast<ScalarType>(280)),
           ScalarType(ONE / static_cast<ScalarType>(14)), ScalarType(ZERO),
           ScalarType(ONE / static_cast<ScalarType>(630)),
           ScalarType(-ONE / static_cast<ScalarType>(1890)),
           ScalarType(ONE / static_cast<ScalarType>(3780)),
           ScalarType(-ONE / static_cast<ScalarType>(252)),
           ScalarType(ONE / static_cast<ScalarType>(18))});
    }
    auto m0_basis = std::array<ScalarType, 3>(
        {signedDistance<ScalarType>(HALF * (pt_0 + pt_1), a_paraboloid),
         signedDistance<ScalarType>(QUARTER * (pt_0 + pt_1) + HALF * cp,
                                    a_paraboloid),
         signedDistance<ScalarType>(cp, a_paraboloid)});
    auto m1x_basis = std::array<ScalarType, 4>(
        {-SIX * (X0Z0 - X0Z2 - X2Z0 + X2Z2 - TWO * B * X2 * Y00 +
                 TWO * B * X0 * Y02 + TWO * B * X2 * Y02 - TWO * B * X0 * Y22),
         TWO *
             (FIVE * X0Z0 + static_cast<ScalarType>(10) * X0Z1p + SIX * X0Z1P +
              static_cast<ScalarType>(7) * X0Z2 +
              static_cast<ScalarType>(30) * A * X02 * X1 -
              static_cast<ScalarType>(11) * X1Z0 - FOUR * X1Z1p -
              static_cast<ScalarType>(11) * X1Z2 +
              static_cast<ScalarType>(7) * X2Z0 +
              static_cast<ScalarType>(10) * X2Z1p + SIX * X2Z1P + FIVE * X2Z2 -
              static_cast<ScalarType>(14) * B * X1 * Y00 + FOUR * B * X2 * Y00 +
              static_cast<ScalarType>(14) * B * X0 * Y01 - FOUR * B * X1 * Y01 +
              static_cast<ScalarType>(10) * B * X2 * Y01 - FOUR * B * X0 * Y02 +
              static_cast<ScalarType>(10) * B * X1 * Y02 - FOUR * B * X2 * Y02 +
              FOUR * B * X0 * Y11 + FOUR * B * X2 * Y11 +
              static_cast<ScalarType>(10) * B * X0 * Y12 - FOUR * B * X1 * Y12 +
              static_cast<ScalarType>(14) * B * X2 * Y12 + FOUR * B * X0 * Y22 -
              static_cast<ScalarType>(14) * B * X1 * Y22),
         TWO *
             (-FIVE * X0Z1p + static_cast<ScalarType>(18) * X0Z1P + X0Z2 +
              SIX * A * X02 * X1 - FIVE * X1Z0 - SIX * X1Z1p - SIX * X1Z1P -
              FIVE * X1Z2 + X2Z0 - FIVE * X2Z1p +
              static_cast<ScalarType>(18) * X2Z1P -
              static_cast<ScalarType>(12) * B * X1 * Y01 + TWO * B * X2 * Y01 +
              TWO * B * X1 * Y02 + static_cast<ScalarType>(12) * B * X0 * Y11 +
              static_cast<ScalarType>(12) * B * X2 * Y11 + TWO * B * X0 * Y12 -
              static_cast<ScalarType>(12) * B * X1 * Y12),
         TWO * (X1Z1p - X1Z1P)});
    auto m1y_basis = std::array<ScalarType, 4>(
        {SIX *
             (-Y0Z0 + Y0Z2 + TWO * A * (X22 * Y0 + X00 * Y2 - X02 * (Y0 + Y2)) +
              Y2Z0 - Y2Z2),
         TWO *
             (FIVE * Y0Z0 + static_cast<ScalarType>(10) * Y0Z1p + SIX * Y0Z1P +
              static_cast<ScalarType>(7) * Y0Z2 +
              static_cast<ScalarType>(30) * B * Y02 * Y1 -
              static_cast<ScalarType>(11) * Y1Z0 - FOUR * Y1Z1p -
              static_cast<ScalarType>(11) * Y1Z2 +
              TWO * A *
                  (-TWO * X02 * Y0 + TWO * X11 * Y0 + FIVE * X12 * Y0 +
                   TWO * X22 * Y0 - static_cast<ScalarType>(7) * X00 * Y1 +
                   FIVE * X02 * Y1 - TWO * X12 * Y1 -
                   static_cast<ScalarType>(7) * X22 * Y1 + TWO * X00 * Y2 -
                   TWO * X02 * Y2 + TWO * X11 * Y2 +
                   static_cast<ScalarType>(7) * X12 * Y2 +
                   X01 * (static_cast<ScalarType>(7) * Y0 - TWO * Y1 +
                          FIVE * Y2)) +
              static_cast<ScalarType>(7) * Y2Z0 +
              static_cast<ScalarType>(10) * Y2Z1p + SIX * Y2Z1P + FIVE * Y2Z2),
         -TWO * (FIVE * Y0Z1p - static_cast<ScalarType>(18) * Y0Z1P - Y0Z2 -
                 SIX * B * Y02 * Y1 + FIVE * Y1Z0 + SIX * Y1Z1p + SIX * Y1Z1P +
                 FIVE * Y1Z2 -
                 TWO * A *
                     (X12 * Y0 - SIX * X01 * Y1 + X02 * Y1 - SIX * X12 * Y1 +
                      X01 * Y2 + SIX * X11 * (Y0 + Y2)) -
                 Y2Z0 + FIVE * Y2Z1p - static_cast<ScalarType>(18) * Y2Z1P),
         TWO * (Y1Z1p - Y1Z1P)});
    auto m1z_basis = std::array<ScalarType, 5>(
        {-(AA * (static_cast<ScalarType>(21) * X0000 +
                 static_cast<ScalarType>(28) * X000 * X2 +
                 static_cast<ScalarType>(30) * X00 * X22 +
                 static_cast<ScalarType>(28) * X0 * X222 +
                 static_cast<ScalarType>(21) * X2222)) -
             static_cast<ScalarType>(21) * BB * Y0000 -
             static_cast<ScalarType>(28) * BB * Y000 * Y2 -
             static_cast<ScalarType>(30) * BB * Y00 * Y22 -
             TWO * AB *
                 (X00 * (static_cast<ScalarType>(21) * Y00 +
                         static_cast<ScalarType>(14) * Y02 + FIVE * Y22) +
                  TWO * X02 *
                      (static_cast<ScalarType>(7) * Y00 +
                       static_cast<ScalarType>(10) * Y02 +
                       static_cast<ScalarType>(7) * Y22) +
                  X22 * (static_cast<ScalarType>(5) * Y00 +
                         static_cast<ScalarType>(14) * Y02 +
                         static_cast<ScalarType>(21) * Y22)) -
             static_cast<ScalarType>(28) * BB * Y0 * Y222 -
             static_cast<ScalarType>(21) * BB * Y2222 +
             static_cast<ScalarType>(40) * Z00 +
             static_cast<ScalarType>(48) * Z02 +
             static_cast<ScalarType>(40) * Z22,
         THREE * AA *
                 (static_cast<ScalarType>(21) * X000 * X1 -
                  static_cast<ScalarType>(7) * X00 * X11 -
                  static_cast<ScalarType>(10) * X02 * X11 +
                  static_cast<ScalarType>(35) * X00 * X12 -
                  static_cast<ScalarType>(7) * X000 * X2 -
                  static_cast<ScalarType>(10) * X00 * X22 +
                  static_cast<ScalarType>(35) * X01 * X22 -
                  static_cast<ScalarType>(7) * X11 * X22 -
                  static_cast<ScalarType>(7) * X0 * X222 +
                  static_cast<ScalarType>(21) * X1 * X222) -
             AB * (static_cast<ScalarType>(7) * X11 * Y00 -
                   static_cast<ScalarType>(35) * X12 * Y00 +
                   static_cast<ScalarType>(10) * X22 * Y00 -
                   static_cast<ScalarType>(63) * X00 * Y01 +
                   static_cast<ScalarType>(20) * X12 * Y01 -
                   static_cast<ScalarType>(35) * X22 * Y01 +
                   static_cast<ScalarType>(21) * X00 * Y02 +
                   static_cast<ScalarType>(10) * X11 * Y02 -
                   static_cast<ScalarType>(70) * X12 * Y02 +
                   static_cast<ScalarType>(21) * X22 * Y02 +
                   static_cast<ScalarType>(7) * X00 * Y11 +
                   static_cast<ScalarType>(7) * X22 * Y11 -
                   static_cast<ScalarType>(35) * X00 * Y12 +
                   static_cast<ScalarType>(28) * X12 * Y12 -
                   static_cast<ScalarType>(63) * X22 * Y12 +
                   X01 * (-static_cast<ScalarType>(63) * Y00 +
                          static_cast<ScalarType>(28) * Y01 -
                          static_cast<ScalarType>(70) * Y02 +
                          static_cast<ScalarType>(20) * Y12 -
                          static_cast<ScalarType>(35) * Y22) +
                   static_cast<ScalarType>(10) * X00 * Y22 +
                   static_cast<ScalarType>(7) * X11 * Y22 -
                   static_cast<ScalarType>(63) * X12 * Y22 +
                   X02 * (static_cast<ScalarType>(21) * Y00 -
                          static_cast<ScalarType>(70) * Y01 +
                          static_cast<ScalarType>(40) * Y02 +
                          static_cast<ScalarType>(10) * Y11 -
                          static_cast<ScalarType>(70) * Y12 +
                          static_cast<ScalarType>(21) * Y22)) -
             THREE *
                 (BB * (static_cast<ScalarType>(7) * Y00 * Y11 +
                        static_cast<ScalarType>(10) * Y02 * Y11 -
                        static_cast<ScalarType>(35) * Y00 * Y12 +
                        static_cast<ScalarType>(7) * Y000 * (-THREE * Y1 + Y2) +
                        static_cast<ScalarType>(10) * Y00 * Y22 -
                        static_cast<ScalarType>(35) * Y01 * Y22 +
                        static_cast<ScalarType>(7) * Y11 * Y22 +
                        static_cast<ScalarType>(7) * Y0 * Y222 -
                        static_cast<ScalarType>(21) * Y1 * Y222) +
                  TWO * (static_cast<ScalarType>(5) * Z00 +
                         static_cast<ScalarType>(10) * Z01p + FOUR * Z02 -
                         TWO * Z1p1p + static_cast<ScalarType>(10) * Z1p2 +
                         FIVE * Z22)),
         -SIX * AA *
                 (static_cast<ScalarType>(46) * X02 * X11 -
                  static_cast<ScalarType>(14) * X0 * X111 + X1111 -
                  static_cast<ScalarType>(14) * X111 * X2 -
                  static_cast<ScalarType>(14) * X01 * X22 +
                  static_cast<ScalarType>(28) * X11 * X22 +
                  X00 * (static_cast<ScalarType>(28) * X11 -
                         static_cast<ScalarType>(14) * X12 + X22)) -
             TWO * AB *
                 (X22 * Y00 + static_cast<ScalarType>(112) * X01 * Y01 -
                  static_cast<ScalarType>(28) * X02 * Y01 -
                  static_cast<ScalarType>(14) * X22 * Y01 -
                  static_cast<ScalarType>(28) * X01 * Y02 + FOUR * X02 * Y02 +
                  static_cast<ScalarType>(28) * X00 * Y11 -
                  static_cast<ScalarType>(42) * X01 * Y11 +
                  static_cast<ScalarType>(46) * X02 * Y11 +
                  static_cast<ScalarType>(28) * X22 * Y11 -
                  TWO * X12 *
                      (static_cast<ScalarType>(7) * Y00 -
                       static_cast<ScalarType>(46) * Y01 +
                       static_cast<ScalarType>(14) * Y02 +
                       static_cast<ScalarType>(21) * Y11 -
                       static_cast<ScalarType>(56) * Y12) -
                  static_cast<ScalarType>(14) * X00 * Y12 +
                  static_cast<ScalarType>(92) * X01 * Y12 -
                  static_cast<ScalarType>(28) * X02 * Y12 + X00 * Y22 -
                  static_cast<ScalarType>(14) * X01 * Y22 +
                  X11 * (static_cast<ScalarType>(28) * Y00 -
                         static_cast<ScalarType>(42) * Y01 +
                         static_cast<ScalarType>(46) * Y02 + SIX * Y11 -
                         static_cast<ScalarType>(42) * Y12 +
                         static_cast<ScalarType>(28) * Y22)) +
             THREE * (-TWO * BB *
                          (static_cast<ScalarType>(46) * Y02 * Y11 -
                           static_cast<ScalarType>(14) * Y0 * Y111 + Y1111 -
                           static_cast<ScalarType>(14) * Y111 * Y2 -
                           static_cast<ScalarType>(14) * Y01 * Y22 +
                           static_cast<ScalarType>(28) * Y11 * Y22 +
                           Y00 * (static_cast<ScalarType>(28) * Y11 -
                                  static_cast<ScalarType>(14) * Y12 + Y22)) +
                      FIVE * Z00 + static_cast<ScalarType>(40) * Z01p -
                      TWO * Z02 + static_cast<ScalarType>(8) * Z1p1p +
                      static_cast<ScalarType>(40) * Z1p2 + FIVE * Z22),
         -TWO * AA *
                 (static_cast<ScalarType>(3) * X02 * X11 -
                  static_cast<ScalarType>(7) * X0 * X111 + THREE * X1111 -
                  static_cast<ScalarType>(7) * X111 * X2) -
             SIX * BB * Y02 * Y11 +
             static_cast<ScalarType>(14) * BB * Y0 * Y111 - SIX * BB * Y1111 -
             TWO * AB *
                 (static_cast<ScalarType>(2) * X12 * Y01 -
                  static_cast<ScalarType>(7) * X01 * Y11 + X02 * Y11 -
                  static_cast<ScalarType>(7) * X12 * Y11 +
                  X11 * (-static_cast<ScalarType>(7) * Y01 + Y02 + SIX * Y11 -
                         static_cast<ScalarType>(7) * Y12) +
                  TWO * X01 * Y12) +
             static_cast<ScalarType>(14) * BB * Y111 * Y2 - FIVE * Z01p + Z02 -
             static_cast<ScalarType>(7) * Z1p1p - FIVE * Z1p2,
         -(AA * X1111) - TWO * AB * X11 * Y11 - BB * Y1111 + Z1p1p});
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
}

// template <class ScalarType>
// inline enable_if_t<!has_embedded_gradient<Volume>::value, Volume>
// computeType3Contribution(const AlignedParaboloidBase<ScalarType>&
// a_paraboloid,
//                          const RationalBezierArcBase<ScalarType>& a_arc) {
//   const auto& pt_0 = a_arc.start_point();
//   const auto& cp = a_arc.control_point();
//   const auto& pt_1 = a_arc.end_point();
//   const auto& weight = a_arc.weight();
//   const ScalarType area_proj_triangle =
//       (pt_0[0] * (pt_1[1] - cp[1]) + pt_1[0] * (cp[1] - pt_0[1]) +
//        cp[0] * (pt_0[1] - pt_1[1])) /
//       static_cast<ScalarType>(2);
//   assert(weight >= static_cast<ScalarType>(0.0));
//   std::array<ScalarType, 3> coeffs;
//   if (weight < static_cast<ScalarType>(0.35))  // We use the exact
//   expressions
//     coeffs = coeffsV3Exact<ScalarType, ScalarType>(weight);
//   else if (weight < static_cast<ScalarType>(
//                         1.7))  // We use the 40th order Taylor series (w ->
//                         1)
//     coeffs = coeffsV3SeriesOne<ScalarType, ScalarType>(weight);
//   else if (weight < static_cast<ScalarType>(
//                         1.0e9))  // We use the series expansion (w ->
//                         infty)
//     coeffs = coeffsV3Exact<ScalarType, ScalarType>(weight);
//   else  // This is within EPSILON of the actual value
//     coeffs = std::array<ScalarType, 3>(
//         {static_cast<ScalarType>(1) / static_cast<ScalarType>(6),
//          static_cast<ScalarType>(2) / static_cast<ScalarType>(3),
//          static_cast<ScalarType>(0)});
//   return area_proj_triangle *
//          (coeffs[0] *
//               signedDistance<ScalarType>(
//                   (pt_0 + pt_1) / static_cast<ScalarType>(2), a_paraboloid)
//                   +
//           coeffs[1] * signedDistance<ScalarType>(
//                           (pt_0 + pt_1) / static_cast<ScalarType>(4) +
//                               cp / static_cast<ScalarType>(2),
//                           a_paraboloid) +
//           coeffs[2] * signedDistance<ScalarType>(cp, a_paraboloid));
// }

template <class ReturnType, class ScalarType>
inline ReturnType computeType3ContributionWithGradient(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const RationalBezierArcWithGradientBase<
        PtWithGradient<typename ReturnType::gradient_type>, ScalarType>&
        a_arc) {
  /* Defining constants and types */
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType THREE = static_cast<ScalarType>(3);
  const ScalarType FOUR = static_cast<ScalarType>(4);
  const ScalarType FIVE = static_cast<ScalarType>(5);
  const ScalarType SIX = static_cast<ScalarType>(6);
  const ScalarType HALF = ONE / TWO;
  const ScalarType QUARTER = ONE / FOUR;

  /* Function */
  using gradient_type = typename ReturnType::gradient_type;
  using scalar_container_type = ScalarWithGradient<gradient_type>;
  auto moments_with_gradient = ReturnType::fromScalarConstant(ZERO);
  std::array<scalar_container_type, 3> pt_0, cp, pt_1;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_0[d] = scalar_container_type(a_arc.start_point().getPt()[d],
                                    a_arc.start_point().getData()[d]);
    cp[d] = scalar_container_type(a_arc.control_point().getPt()[d],
                                  a_arc.control_point().getData()[d]);
    pt_1[d] = scalar_container_type(a_arc.end_point().getPt()[d],
                                    a_arc.end_point().getData()[d]);
  }
  auto A_grad = gradient_type(ZERO), B_grad = gradient_type(ZERO);
  A_grad.setGradA(ONE);
  B_grad.setGradB(ONE);
  const auto A = scalar_container_type(a_paraboloid.a(), A_grad),
             B = scalar_container_type(a_paraboloid.b(), B_grad);
  const auto X0 = pt_0[0], X1 = cp[0], X2 = pt_1[0];
  const auto Y0 = pt_0[1], Y1 = cp[1], Y2 = pt_1[1];
  const auto Z0 = pt_0[2], Z1p = cp[2], Z2 = pt_1[2];
  const auto area_proj_triangle =
      HALF * (X0 * (Y2 - Y1) + X2 * (Y1 - Y0) + X1 * (Y0 - Y2));
  const auto m0_basis = std::array<scalar_container_type, 3>(
      {area_proj_triangle *
           (QUARTER * A * (X0 + X2) * (X0 + X2) +
            QUARTER * B * (Y0 + Y2) * (Y0 + Y2) + HALF * (Z0 + Z2)),
       area_proj_triangle * (A * (QUARTER * (X0 + X2) + HALF * X1) *
                                 (QUARTER * (X0 + X2) + HALF * X1) +
                             B * (QUARTER * (Y0 + Y2) + HALF * Y1) *
                                 (QUARTER * (Y0 + Y2) + HALF * Y1) +
                             (QUARTER * (Z0 + Z2) + HALF * Z1p)),
       area_proj_triangle * (A * X1 * X1 + B * Y1 * Y1 + Z1p)});
  const auto weight =
      scalar_container_type(a_arc.weight(), a_arc.weight_gradient());
  assert(weight.value() >= ZERO);
  if constexpr (std::is_same<typename ReturnType::moment_type, Volume>::value) {
    std::array<scalar_container_type, 3> coeffs;
    if (weight < scalar_container_type(static_cast<ScalarType>(
                     0.35)))  // We use the exact expressions
      coeffs = coeffsV3Exact<scalar_container_type, ScalarType>(weight);
    else if (weight <
             scalar_container_type(static_cast<ScalarType>(
                 1.7)))  // We use the 40th order Taylor series (w -> 1)
      coeffs = coeffsV3SeriesOne<scalar_container_type, ScalarType>(weight);
    else if (weight < scalar_container_type(static_cast<ScalarType>(
                          1.0e9)))  // We use the series expansion (w -> infty)
      coeffs = coeffsV3Exact<scalar_container_type, ScalarType>(weight);
    else  // This is within EPSILON of the actual value
      coeffs = std::array<scalar_container_type, 3>(
          {scalar_container_type(ONE / SIX), scalar_container_type(TWO / THREE),
           scalar_container_type(ZERO)});
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
    std::array<scalar_container_type, 12> coeffs;
    if (weight < scalar_container_type(static_cast<ScalarType>(
                     0.35)))  // We use the exact expressions
    {
      coeffs = coeffsV3andC3Exact<scalar_container_type, ScalarType>(weight);
    } else if (weight < scalar_container_type(static_cast<ScalarType>(
                            1.7)))  // We use Taylor series (w -> 1)
    {
      coeffs =
          coeffsV3andC3SeriesOne<scalar_container_type, ScalarType>(weight);
    } else if (weight < scalar_container_type(static_cast<ScalarType>(
                            1.0e9)))  // We use the exact expressions
    {
      coeffs = coeffsV3andC3Exact<scalar_container_type, ScalarType>(weight);
    }
    // else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
    //   coeffs = coeffsV3andC3SeriesInfinity(weight);
    else  // This is within EPSILON of the actual value
    {
      coeffs = std::array<scalar_container_type, 12>(
          {scalar_container_type(ONE / SIX), scalar_container_type(TWO / THREE),
           scalar_container_type(-ONE / static_cast<ScalarType>(140)),
           scalar_container_type(ONE / static_cast<ScalarType>(420)),
           scalar_container_type(-ONE / static_cast<ScalarType>(280)),
           scalar_container_type(ONE / static_cast<ScalarType>(14)),
           scalar_container_type(ZERO),
           scalar_container_type(ONE / static_cast<ScalarType>(630)),
           scalar_container_type(-ONE / static_cast<ScalarType>(1890)),
           scalar_container_type(ONE / static_cast<ScalarType>(3780)),
           scalar_container_type(-ONE / static_cast<ScalarType>(252)),
           scalar_container_type(ONE / static_cast<ScalarType>(18))});
    }
    auto m1x_basis = std::array<scalar_container_type, 4>(
        {-SIX * (X0Z0 - X0Z2 - X2Z0 + X2Z2 - TWO * B * X2 * Y00 +
                 TWO * B * X0 * Y02 + TWO * B * X2 * Y02 - TWO * B * X0 * Y22),
         TWO *
             (FIVE * X0Z0 + static_cast<ScalarType>(10) * X0Z1p + SIX * X0Z1P +
              static_cast<ScalarType>(7) * X0Z2 +
              static_cast<ScalarType>(30) * A * X02 * X1 -
              static_cast<ScalarType>(11) * X1Z0 - FOUR * X1Z1p -
              static_cast<ScalarType>(11) * X1Z2 +
              static_cast<ScalarType>(7) * X2Z0 +
              static_cast<ScalarType>(10) * X2Z1p + SIX * X2Z1P + FIVE * X2Z2 -
              static_cast<ScalarType>(14) * B * X1 * Y00 + FOUR * B * X2 * Y00 +
              static_cast<ScalarType>(14) * B * X0 * Y01 - FOUR * B * X1 * Y01 +
              static_cast<ScalarType>(10) * B * X2 * Y01 - FOUR * B * X0 * Y02 +
              static_cast<ScalarType>(10) * B * X1 * Y02 - FOUR * B * X2 * Y02 +
              FOUR * B * X0 * Y11 + FOUR * B * X2 * Y11 +
              static_cast<ScalarType>(10) * B * X0 * Y12 - FOUR * B * X1 * Y12 +
              static_cast<ScalarType>(14) * B * X2 * Y12 + FOUR * B * X0 * Y22 -
              static_cast<ScalarType>(14) * B * X1 * Y22),
         TWO *
             (-FIVE * X0Z1p + static_cast<ScalarType>(18) * X0Z1P + X0Z2 +
              SIX * A * X02 * X1 - FIVE * X1Z0 - SIX * X1Z1p - SIX * X1Z1P -
              FIVE * X1Z2 + X2Z0 - FIVE * X2Z1p +
              static_cast<ScalarType>(18) * X2Z1P -
              static_cast<ScalarType>(12) * B * X1 * Y01 + TWO * B * X2 * Y01 +
              TWO * B * X1 * Y02 + static_cast<ScalarType>(12) * B * X0 * Y11 +
              static_cast<ScalarType>(12) * B * X2 * Y11 + TWO * B * X0 * Y12 -
              static_cast<ScalarType>(12) * B * X1 * Y12),
         TWO * (X1Z1p - X1Z1P)});
    auto m1y_basis = std::array<scalar_container_type, 4>(
        {SIX *
             (-Y0Z0 + Y0Z2 + TWO * A * (X22 * Y0 + X00 * Y2 - X02 * (Y0 + Y2)) +
              Y2Z0 - Y2Z2),
         TWO *
             (FIVE * Y0Z0 + static_cast<ScalarType>(10) * Y0Z1p + SIX * Y0Z1P +
              static_cast<ScalarType>(7) * Y0Z2 +
              static_cast<ScalarType>(30) * B * Y02 * Y1 -
              static_cast<ScalarType>(11) * Y1Z0 - FOUR * Y1Z1p -
              static_cast<ScalarType>(11) * Y1Z2 +
              TWO * A *
                  (-TWO * X02 * Y0 + TWO * X11 * Y0 + FIVE * X12 * Y0 +
                   TWO * X22 * Y0 - static_cast<ScalarType>(7) * X00 * Y1 +
                   FIVE * X02 * Y1 - TWO * X12 * Y1 -
                   static_cast<ScalarType>(7) * X22 * Y1 + TWO * X00 * Y2 -
                   TWO * X02 * Y2 + TWO * X11 * Y2 +
                   static_cast<ScalarType>(7) * X12 * Y2 +
                   X01 * (static_cast<ScalarType>(7) * Y0 - TWO * Y1 +
                          FIVE * Y2)) +
              static_cast<ScalarType>(7) * Y2Z0 +
              static_cast<ScalarType>(10) * Y2Z1p + SIX * Y2Z1P + FIVE * Y2Z2),
         -TWO * (FIVE * Y0Z1p - static_cast<ScalarType>(18) * Y0Z1P - Y0Z2 -
                 SIX * B * Y02 * Y1 + FIVE * Y1Z0 + SIX * Y1Z1p + SIX * Y1Z1P +
                 FIVE * Y1Z2 -
                 TWO * A *
                     (X12 * Y0 - SIX * X01 * Y1 + X02 * Y1 - SIX * X12 * Y1 +
                      X01 * Y2 + SIX * X11 * (Y0 + Y2)) -
                 Y2Z0 + FIVE * Y2Z1p - static_cast<ScalarType>(18) * Y2Z1P),
         TWO * (Y1Z1p - Y1Z1P)});
    auto m1z_basis = std::array<scalar_container_type, 5>(
        {-(AA * (static_cast<ScalarType>(21) * X0000 +
                 static_cast<ScalarType>(28) * X000 * X2 +
                 static_cast<ScalarType>(30) * X00 * X22 +
                 static_cast<ScalarType>(28) * X0 * X222 +
                 static_cast<ScalarType>(21) * X2222)) -
             static_cast<ScalarType>(21) * BB * Y0000 -
             static_cast<ScalarType>(28) * BB * Y000 * Y2 -
             static_cast<ScalarType>(30) * BB * Y00 * Y22 -
             TWO * AB *
                 (X00 * (static_cast<ScalarType>(21) * Y00 +
                         static_cast<ScalarType>(14) * Y02 + FIVE * Y22) +
                  TWO * X02 *
                      (static_cast<ScalarType>(7) * Y00 +
                       static_cast<ScalarType>(10) * Y02 +
                       static_cast<ScalarType>(7) * Y22) +
                  X22 * (static_cast<ScalarType>(5) * Y00 +
                         static_cast<ScalarType>(14) * Y02 +
                         static_cast<ScalarType>(21) * Y22)) -
             static_cast<ScalarType>(28) * BB * Y0 * Y222 -
             static_cast<ScalarType>(21) * BB * Y2222 +
             static_cast<ScalarType>(40) * Z00 +
             static_cast<ScalarType>(48) * Z02 +
             static_cast<ScalarType>(40) * Z22,
         THREE * AA *
                 (static_cast<ScalarType>(21) * X000 * X1 -
                  static_cast<ScalarType>(7) * X00 * X11 -
                  static_cast<ScalarType>(10) * X02 * X11 +
                  static_cast<ScalarType>(35) * X00 * X12 -
                  static_cast<ScalarType>(7) * X000 * X2 -
                  static_cast<ScalarType>(10) * X00 * X22 +
                  static_cast<ScalarType>(35) * X01 * X22 -
                  static_cast<ScalarType>(7) * X11 * X22 -
                  static_cast<ScalarType>(7) * X0 * X222 +
                  static_cast<ScalarType>(21) * X1 * X222) -
             AB * (static_cast<ScalarType>(7) * X11 * Y00 -
                   static_cast<ScalarType>(35) * X12 * Y00 +
                   static_cast<ScalarType>(10) * X22 * Y00 -
                   static_cast<ScalarType>(63) * X00 * Y01 +
                   static_cast<ScalarType>(20) * X12 * Y01 -
                   static_cast<ScalarType>(35) * X22 * Y01 +
                   static_cast<ScalarType>(21) * X00 * Y02 +
                   static_cast<ScalarType>(10) * X11 * Y02 -
                   static_cast<ScalarType>(70) * X12 * Y02 +
                   static_cast<ScalarType>(21) * X22 * Y02 +
                   static_cast<ScalarType>(7) * X00 * Y11 +
                   static_cast<ScalarType>(7) * X22 * Y11 -
                   static_cast<ScalarType>(35) * X00 * Y12 +
                   static_cast<ScalarType>(28) * X12 * Y12 -
                   static_cast<ScalarType>(63) * X22 * Y12 +
                   X01 * (-static_cast<ScalarType>(63) * Y00 +
                          static_cast<ScalarType>(28) * Y01 -
                          static_cast<ScalarType>(70) * Y02 +
                          static_cast<ScalarType>(20) * Y12 -
                          static_cast<ScalarType>(35) * Y22) +
                   static_cast<ScalarType>(10) * X00 * Y22 +
                   static_cast<ScalarType>(7) * X11 * Y22 -
                   static_cast<ScalarType>(63) * X12 * Y22 +
                   X02 * (static_cast<ScalarType>(21) * Y00 -
                          static_cast<ScalarType>(70) * Y01 +
                          static_cast<ScalarType>(40) * Y02 +
                          static_cast<ScalarType>(10) * Y11 -
                          static_cast<ScalarType>(70) * Y12 +
                          static_cast<ScalarType>(21) * Y22)) -
             THREE *
                 (BB * (static_cast<ScalarType>(7) * Y00 * Y11 +
                        static_cast<ScalarType>(10) * Y02 * Y11 -
                        static_cast<ScalarType>(35) * Y00 * Y12 +
                        static_cast<ScalarType>(7) * Y000 * (-THREE * Y1 + Y2) +
                        static_cast<ScalarType>(10) * Y00 * Y22 -
                        static_cast<ScalarType>(35) * Y01 * Y22 +
                        static_cast<ScalarType>(7) * Y11 * Y22 +
                        static_cast<ScalarType>(7) * Y0 * Y222 -
                        static_cast<ScalarType>(21) * Y1 * Y222) +
                  TWO * (static_cast<ScalarType>(5) * Z00 +
                         static_cast<ScalarType>(10) * Z01p + FOUR * Z02 -
                         TWO * Z1p1p + static_cast<ScalarType>(10) * Z1p2 +
                         FIVE * Z22)),
         -SIX * AA *
                 (static_cast<ScalarType>(46) * X02 * X11 -
                  static_cast<ScalarType>(14) * X0 * X111 + X1111 -
                  static_cast<ScalarType>(14) * X111 * X2 -
                  static_cast<ScalarType>(14) * X01 * X22 +
                  static_cast<ScalarType>(28) * X11 * X22 +
                  X00 * (static_cast<ScalarType>(28) * X11 -
                         static_cast<ScalarType>(14) * X12 + X22)) -
             TWO * AB *
                 (X22 * Y00 + static_cast<ScalarType>(112) * X01 * Y01 -
                  static_cast<ScalarType>(28) * X02 * Y01 -
                  static_cast<ScalarType>(14) * X22 * Y01 -
                  static_cast<ScalarType>(28) * X01 * Y02 + FOUR * X02 * Y02 +
                  static_cast<ScalarType>(28) * X00 * Y11 -
                  static_cast<ScalarType>(42) * X01 * Y11 +
                  static_cast<ScalarType>(46) * X02 * Y11 +
                  static_cast<ScalarType>(28) * X22 * Y11 -
                  TWO * X12 *
                      (static_cast<ScalarType>(7) * Y00 -
                       static_cast<ScalarType>(46) * Y01 +
                       static_cast<ScalarType>(14) * Y02 +
                       static_cast<ScalarType>(21) * Y11 -
                       static_cast<ScalarType>(56) * Y12) -
                  static_cast<ScalarType>(14) * X00 * Y12 +
                  static_cast<ScalarType>(92) * X01 * Y12 -
                  static_cast<ScalarType>(28) * X02 * Y12 + X00 * Y22 -
                  static_cast<ScalarType>(14) * X01 * Y22 +
                  X11 * (static_cast<ScalarType>(28) * Y00 -
                         static_cast<ScalarType>(42) * Y01 +
                         static_cast<ScalarType>(46) * Y02 + SIX * Y11 -
                         static_cast<ScalarType>(42) * Y12 +
                         static_cast<ScalarType>(28) * Y22)) +
             THREE * (-TWO * BB *
                          (static_cast<ScalarType>(46) * Y02 * Y11 -
                           static_cast<ScalarType>(14) * Y0 * Y111 + Y1111 -
                           static_cast<ScalarType>(14) * Y111 * Y2 -
                           static_cast<ScalarType>(14) * Y01 * Y22 +
                           static_cast<ScalarType>(28) * Y11 * Y22 +
                           Y00 * (static_cast<ScalarType>(28) * Y11 -
                                  static_cast<ScalarType>(14) * Y12 + Y22)) +
                      FIVE * Z00 + static_cast<ScalarType>(40) * Z01p -
                      TWO * Z02 + static_cast<ScalarType>(8) * Z1p1p +
                      static_cast<ScalarType>(40) * Z1p2 + FIVE * Z22),
         -TWO * AA *
                 (static_cast<ScalarType>(3) * X02 * X11 -
                  static_cast<ScalarType>(7) * X0 * X111 + THREE * X1111 -
                  static_cast<ScalarType>(7) * X111 * X2) -
             SIX * BB * Y02 * Y11 +
             static_cast<ScalarType>(14) * BB * Y0 * Y111 - SIX * BB * Y1111 -
             TWO * AB *
                 (static_cast<ScalarType>(2) * X12 * Y01 -
                  static_cast<ScalarType>(7) * X01 * Y11 + X02 * Y11 -
                  static_cast<ScalarType>(7) * X12 * Y11 +
                  X11 * (-static_cast<ScalarType>(7) * Y01 + Y02 + SIX * Y11 -
                         static_cast<ScalarType>(7) * Y12) +
                  TWO * X01 * Y12) +
             static_cast<ScalarType>(14) * BB * Y111 * Y2 - FIVE * Z01p + Z02 -
             static_cast<ScalarType>(7) * Z1p1p - FIVE * Z1p2,
         -(AA * X1111) - TWO * AB * X11 * Y11 - BB * Y1111 + Z1p1p});
    // Update Volume
    auto m0 = scalar_container_type(ZERO);
    for (size_t i = 0; i < 3; ++i) {
      m0 += coeffs[i] * m0_basis[i];
    }
    // Update Centroid
    auto m1x = scalar_container_type(ZERO), m1y = scalar_container_type(ZERO),
         m1z = scalar_container_type(ZERO);
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
    PtBase<ScalarType>& centroid = moments_with_gradient.centroid().getPt();
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

// template <class ScalarType>
// inline enable_if_t<!has_embedded_gradient<VolumeMoments>::value,
// VolumeMoments> computeType3Contribution(const
// AlignedParaboloidBase<ScalarType>& a_paraboloid,
//                          const RationalBezierArcBase<ScalarType>& a_arc) {
//   /* Defining constants and types */
//   const ScalarType ZERO = static_cast<ScalarType>(0);
//   const ScalarType ONE = static_cast<ScalarType>(1);
//   const ScalarType TWO = static_cast<ScalarType>(2);
//   const ScalarType THREE = static_cast<ScalarType>(3);
//   const ScalarType FOUR = static_cast<ScalarType>(4);
//   const ScalarType FIVE = static_cast<ScalarType>(5);
//   const ScalarType SIX = static_cast<ScalarType>(6);
//   const ScalarType HALF = ONE / TWO;
//   const ScalarType QUARTER = ONE / FOUR;

//   /* Function */
//   auto moments = VolumeMoments::fromScalarConstant(ZERO);
//   const auto& pt_0 = a_arc.start_point();
//   const auto& cp = a_arc.control_point();
//   const auto& pt_1 = a_arc.end_point();
//   const auto weight = a_arc.weight();
//   const auto A = a_paraboloid.a(), B = a_paraboloid.b();
//   const auto X0 = pt_0[0], X1 = cp[0], X2 = pt_1[0];
//   const auto Y0 = pt_0[1], Y1 = cp[1], Y2 = pt_1[1];
//   const auto Z0 = pt_0[2], Z1p = cp[2], Z2 = pt_1[2];
//   const ScalarType Z1P = -A * X1 * X1 - B * Y1 * Y1;
//   const ScalarType AA = A * A, BB = B * B, AB = A * B;
//   const ScalarType X00 = X0 * X0, X11 = X1 * X1, X22 = X2 * X2;
//   const ScalarType X000 = X00 * X0, X111 = X11 * X1, X222 = X22 * X2;
//   const ScalarType X0000 = X00 * X00, X1111 = X11 * X11, X2222 = X22 * X22;
//   const ScalarType Y00 = Y0 * Y0, Y11 = Y1 * Y1, Y22 = Y2 * Y2;
//   const ScalarType Y000 = Y00 * Y0, Y111 = Y11 * Y1, Y222 = Y22 * Y2;
//   const ScalarType Y0000 = Y00 * Y00, Y1111 = Y11 * Y11, Y2222 = Y22 * Y22;
//   const ScalarType Z00 = Z0 * Z0, Z22 = Z2 * Z2;
//   const ScalarType Z1p1p = Z1p * Z1p, Z1P1P = Z1P * Z1P;
//   const ScalarType X02 = X0 * X2, X12 = X1 * X2, X01 = X0 * X1;
//   const ScalarType Y02 = Y0 * Y2, Y12 = Y1 * Y2, Y01 = Y0 * Y1;
//   const ScalarType Z02 = Z0 * Z2, Z01p = Z0 * Z1p, Z01P = Z0 * Z1P;
//   const ScalarType Z1p1P = Z1p * Z1P, Z1p2 = Z1p * Z2, Z1P2 = Z1P * Z2;
//   const ScalarType X0Z0 = X0 * Z0, X0Z1p = X0 * Z1p, X0Z1P = X0 * Z1P,
//                    X0Z2 = X0 * Z2;
//   const ScalarType X1Z0 = X1 * Z0, X1Z1p = X1 * Z1p, X1Z1P = X1 * Z1P,
//                    X1Z2 = X1 * Z2;
//   const ScalarType X2Z0 = X2 * Z0, X2Z1p = X2 * Z1p, X2Z1P = X2 * Z1P,
//                    X2Z2 = X2 * Z2;
//   const ScalarType Y0Z0 = Y0 * Z0, Y0Z1p = Y0 * Z1p, Y0Z1P = Y0 * Z1P,
//                    Y0Z2 = Y0 * Z2;
//   const ScalarType Y1Z0 = Y1 * Z0, Y1Z1p = Y1 * Z1p, Y1Z1P = Y1 * Z1P,
//                    Y1Z2 = Y1 * Z2;
//   const ScalarType Y2Z0 = Y2 * Z0, Y2Z1p = Y2 * Z1p, Y2Z1P = Y2 * Z1P,
//                    Y2Z2 = Y2 * Z2;
//   const ScalarType area_proj_triangle =
//       HALF * (X0 * (Y2 - Y1) + X2 * (Y1 - Y0) + X1 * (Y0 - Y2));
//   assert(weight >= ZERO);
//   // Compute coefficients (functions of the weight)
//   std::array<ScalarType, 12> coeffs;
//   if (weight < static_cast<ScalarType>(0.35))  // We use the exact
//   expressions
//   {
//     coeffs = coeffsV3andC3Exact<ScalarType, ScalarType>(weight);
//   } else if (weight < static_cast<ScalarType>(
//                           1.7))  // We use the 40th order Taylor series (w
//                           -> 1)
//   {
//     coeffs = coeffsV3andC3SeriesOne<ScalarType, ScalarType>(weight);
//   } else if (weight <
//              static_cast<ScalarType>(1.0e9))  // We use the exact
//              expressions
//   {
//     coeffs = coeffsV3andC3Exact<ScalarType, ScalarType>(weight);
//   }
//   // else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
//   //   coeffs = coeffsV3andC3SeriesInfinity(weight);
//   else  // This is within EPSILON of the actual value
//   {
//     coeffs = std::array<ScalarType, 12>(
//         {ScalarType(ONE / SIX), ScalarType(TWO / THREE),
//          ScalarType(-ONE / static_cast<ScalarType>(140)),
//          ScalarType(ONE / static_cast<ScalarType>(420)),
//          ScalarType(-ONE / static_cast<ScalarType>(280)),
//          ScalarType(ONE / static_cast<ScalarType>(14)), ScalarType(ZERO),
//          ScalarType(ONE / static_cast<ScalarType>(630)),
//          ScalarType(-ONE / static_cast<ScalarType>(1890)),
//          ScalarType(ONE / static_cast<ScalarType>(3780)),
//          ScalarType(-ONE / static_cast<ScalarType>(252)),
//          ScalarType(ONE / static_cast<ScalarType>(18))});
//   }
//   auto m0_basis = std::array<ScalarType, 3>(
//       {signedDistance<ScalarType>(HALF * (pt_0 + pt_1), a_paraboloid),
//        signedDistance<ScalarType>(QUARTER * (pt_0 + pt_1) + HALF * cp,
//                                   a_paraboloid),
//        signedDistance<ScalarType>(cp, a_paraboloid)});
//   auto m1x_basis = std::array<ScalarType, 4>(
//       {-SIX * (X0Z0 - X0Z2 - X2Z0 + X2Z2 - TWO * B * X2 * Y00 +
//                TWO * B * X0 * Y02 + TWO * B * X2 * Y02 - TWO * B * X0 *
//                Y22),
//        TWO * (FIVE * X0Z0 + static_cast<ScalarType>(10) * X0Z1p + SIX *
//        X0Z1P
//        +
//               static_cast<ScalarType>(7) * X0Z2 +
//               static_cast<ScalarType>(30) * A * X02 * X1 -
//               static_cast<ScalarType>(11) * X1Z0 - FOUR * X1Z1p -
//               static_cast<ScalarType>(11) * X1Z2 +
//               static_cast<ScalarType>(7) * X2Z0 +
//               static_cast<ScalarType>(10) * X2Z1p + SIX * X2Z1P + FIVE *
//               X2Z2
//               - static_cast<ScalarType>(14) * B * X1 * Y00 + FOUR * B * X2
//               * Y00 + static_cast<ScalarType>(14) * B * X0 * Y01 - FOUR * B
//               * X1 * Y01 + static_cast<ScalarType>(10) * B * X2 * Y01 -
//               FOUR * B * X0 * Y02 + static_cast<ScalarType>(10) * B * X1 *
//               Y02 - FOUR * B * X2 * Y02 + FOUR * B * X0 * Y11 + FOUR * B *
//               X2 * Y11
//               + static_cast<ScalarType>(10) * B * X0 * Y12 - FOUR * B * X1
//               * Y12 + static_cast<ScalarType>(14) * B * X2 * Y12 + FOUR * B
//               * X0 * Y22 - static_cast<ScalarType>(14) * B * X1 * Y22),
//        TWO * (-FIVE * X0Z1p + static_cast<ScalarType>(18) * X0Z1P + X0Z2 +
//               SIX * A * X02 * X1 - FIVE * X1Z0 - SIX * X1Z1p - SIX * X1Z1P
//               - FIVE * X1Z2 + X2Z0 - FIVE * X2Z1p +
//               static_cast<ScalarType>(18) * X2Z1P -
//               static_cast<ScalarType>(12) * B * X1 * Y01 + TWO * B * X2 *
//               Y01
//               + TWO * B * X1 * Y02 + static_cast<ScalarType>(12) * B * X0 *
//               Y11 + static_cast<ScalarType>(12) * B * X2 * Y11 + TWO * B *
//               X0
//               * Y12 - static_cast<ScalarType>(12) * B * X1 * Y12),
//        TWO * (X1Z1p - X1Z1P)});
//   auto m1y_basis = std::array<ScalarType, 4>(
//       {SIX * (-Y0Z0 + Y0Z2 + TWO * A * (X22 * Y0 + X00 * Y2 - X02 * (Y0 +
//       Y2)) +
//               Y2Z0 - Y2Z2),
//        TWO * (FIVE * Y0Z0 + static_cast<ScalarType>(10) * Y0Z1p + SIX *
//        Y0Z1P
//        +
//               static_cast<ScalarType>(7) * Y0Z2 +
//               static_cast<ScalarType>(30) * B * Y02 * Y1 -
//               static_cast<ScalarType>(11) * Y1Z0 - FOUR * Y1Z1p -
//               static_cast<ScalarType>(11) * Y1Z2 +
//               TWO * A *
//                   (-TWO * X02 * Y0 + TWO * X11 * Y0 + FIVE * X12 * Y0 +
//                    TWO * X22 * Y0 - static_cast<ScalarType>(7) * X00 * Y1 +
//                    FIVE * X02 * Y1 - TWO * X12 * Y1 -
//                    static_cast<ScalarType>(7) * X22 * Y1 + TWO * X00 * Y2 -
//                    TWO * X02 * Y2 + TWO * X11 * Y2 +
//                    static_cast<ScalarType>(7) * X12 * Y2 +
//                    X01 * (static_cast<ScalarType>(7) * Y0 - TWO * Y1 +
//                           FIVE * Y2)) +
//               static_cast<ScalarType>(7) * Y2Z0 +
//               static_cast<ScalarType>(10) * Y2Z1p + SIX * Y2Z1P + FIVE *
//               Y2Z2),
//        -TWO * (FIVE * Y0Z1p - static_cast<ScalarType>(18) * Y0Z1P - Y0Z2 -
//                SIX * B * Y02 * Y1 + FIVE * Y1Z0 + SIX * Y1Z1p + SIX * Y1Z1P
//                + FIVE * Y1Z2 - TWO * A *
//                    (X12 * Y0 - SIX * X01 * Y1 + X02 * Y1 - SIX * X12 * Y1 +
//                     X01 * Y2 + SIX * X11 * (Y0 + Y2)) -
//                Y2Z0 + FIVE * Y2Z1p - static_cast<ScalarType>(18) * Y2Z1P),
//        TWO * (Y1Z1p - Y1Z1P)});
//   auto m1z_basis = std::array<ScalarType, 5>(
//       {-(AA * (static_cast<ScalarType>(21) * X0000 +
//                static_cast<ScalarType>(28) * X000 * X2 +
//                static_cast<ScalarType>(30) * X00 * X22 +
//                static_cast<ScalarType>(28) * X0 * X222 +
//                static_cast<ScalarType>(21) * X2222)) -
//            static_cast<ScalarType>(21) * BB * Y0000 -
//            static_cast<ScalarType>(28) * BB * Y000 * Y2 -
//            static_cast<ScalarType>(30) * BB * Y00 * Y22 -
//            TWO * AB *
//                (X00 * (static_cast<ScalarType>(21) * Y00 +
//                        static_cast<ScalarType>(14) * Y02 + FIVE * Y22) +
//                 TWO * X02 *
//                     (static_cast<ScalarType>(7) * Y00 +
//                      static_cast<ScalarType>(10) * Y02 +
//                      static_cast<ScalarType>(7) * Y22) +
//                 X22 * (static_cast<ScalarType>(5) * Y00 +
//                        static_cast<ScalarType>(14) * Y02 +
//                        static_cast<ScalarType>(21) * Y22)) -
//            static_cast<ScalarType>(28) * BB * Y0 * Y222 -
//            static_cast<ScalarType>(21) * BB * Y2222 +
//            static_cast<ScalarType>(40) * Z00 +
//            static_cast<ScalarType>(48) * Z02 +
//            static_cast<ScalarType>(40) * Z22,
//        THREE * AA *
//                (static_cast<ScalarType>(21) * X000 * X1 -
//                 static_cast<ScalarType>(7) * X00 * X11 -
//                 static_cast<ScalarType>(10) * X02 * X11 +
//                 static_cast<ScalarType>(35) * X00 * X12 -
//                 static_cast<ScalarType>(7) * X000 * X2 -
//                 static_cast<ScalarType>(10) * X00 * X22 +
//                 static_cast<ScalarType>(35) * X01 * X22 -
//                 static_cast<ScalarType>(7) * X11 * X22 -
//                 static_cast<ScalarType>(7) * X0 * X222 +
//                 static_cast<ScalarType>(21) * X1 * X222) -
//            AB * (static_cast<ScalarType>(7) * X11 * Y00 -
//                  static_cast<ScalarType>(35) * X12 * Y00 +
//                  static_cast<ScalarType>(10) * X22 * Y00 -
//                  static_cast<ScalarType>(63) * X00 * Y01 +
//                  static_cast<ScalarType>(20) * X12 * Y01 -
//                  static_cast<ScalarType>(35) * X22 * Y01 +
//                  static_cast<ScalarType>(21) * X00 * Y02 +
//                  static_cast<ScalarType>(10) * X11 * Y02 -
//                  static_cast<ScalarType>(70) * X12 * Y02 +
//                  static_cast<ScalarType>(21) * X22 * Y02 +
//                  static_cast<ScalarType>(7) * X00 * Y11 +
//                  static_cast<ScalarType>(7) * X22 * Y11 -
//                  static_cast<ScalarType>(35) * X00 * Y12 +
//                  static_cast<ScalarType>(28) * X12 * Y12 -
//                  static_cast<ScalarType>(63) * X22 * Y12 +
//                  X01 * (-static_cast<ScalarType>(63) * Y00 +
//                         static_cast<ScalarType>(28) * Y01 -
//                         static_cast<ScalarType>(70) * Y02 +
//                         static_cast<ScalarType>(20) * Y12 -
//                         static_cast<ScalarType>(35) * Y22) +
//                  static_cast<ScalarType>(10) * X00 * Y22 +
//                  static_cast<ScalarType>(7) * X11 * Y22 -
//                  static_cast<ScalarType>(63) * X12 * Y22 +
//                  X02 * (static_cast<ScalarType>(21) * Y00 -
//                         static_cast<ScalarType>(70) * Y01 +
//                         static_cast<ScalarType>(40) * Y02 +
//                         static_cast<ScalarType>(10) * Y11 -
//                         static_cast<ScalarType>(70) * Y12 +
//                         static_cast<ScalarType>(21) * Y22)) -
//            THREE *
//                (BB * (static_cast<ScalarType>(7) * Y00 * Y11 +
//                       static_cast<ScalarType>(10) * Y02 * Y11 -
//                       static_cast<ScalarType>(35) * Y00 * Y12 +
//                       static_cast<ScalarType>(7) * Y000 * (-THREE * Y1 +
//                       Y2)
//                       + static_cast<ScalarType>(10) * Y00 * Y22 -
//                       static_cast<ScalarType>(35) * Y01 * Y22 +
//                       static_cast<ScalarType>(7) * Y11 * Y22 +
//                       static_cast<ScalarType>(7) * Y0 * Y222 -
//                       static_cast<ScalarType>(21) * Y1 * Y222) +
//                 TWO * (static_cast<ScalarType>(5) * Z00 +
//                        static_cast<ScalarType>(10) * Z01p + FOUR * Z02 -
//                        TWO * Z1p1p + static_cast<ScalarType>(10) * Z1p2 +
//                        FIVE * Z22)),
//        -SIX * AA *
//                (static_cast<ScalarType>(46) * X02 * X11 -
//                 static_cast<ScalarType>(14) * X0 * X111 + X1111 -
//                 static_cast<ScalarType>(14) * X111 * X2 -
//                 static_cast<ScalarType>(14) * X01 * X22 +
//                 static_cast<ScalarType>(28) * X11 * X22 +
//                 X00 * (static_cast<ScalarType>(28) * X11 -
//                        static_cast<ScalarType>(14) * X12 + X22)) -
//            TWO * AB *
//                (X22 * Y00 + static_cast<ScalarType>(112) * X01 * Y01 -
//                 static_cast<ScalarType>(28) * X02 * Y01 -
//                 static_cast<ScalarType>(14) * X22 * Y01 -
//                 static_cast<ScalarType>(28) * X01 * Y02 + FOUR * X02 * Y02
//                 + static_cast<ScalarType>(28) * X00 * Y11 -
//                 static_cast<ScalarType>(42) * X01 * Y11 +
//                 static_cast<ScalarType>(46) * X02 * Y11 +
//                 static_cast<ScalarType>(28) * X22 * Y11 -
//                 TWO * X12 *
//                     (static_cast<ScalarType>(7) * Y00 -
//                      static_cast<ScalarType>(46) * Y01 +
//                      static_cast<ScalarType>(14) * Y02 +
//                      static_cast<ScalarType>(21) * Y11 -
//                      static_cast<ScalarType>(56) * Y12) -
//                 static_cast<ScalarType>(14) * X00 * Y12 +
//                 static_cast<ScalarType>(92) * X01 * Y12 -
//                 static_cast<ScalarType>(28) * X02 * Y12 + X00 * Y22 -
//                 static_cast<ScalarType>(14) * X01 * Y22 +
//                 X11 * (static_cast<ScalarType>(28) * Y00 -
//                        static_cast<ScalarType>(42) * Y01 +
//                        static_cast<ScalarType>(46) * Y02 + SIX * Y11 -
//                        static_cast<ScalarType>(42) * Y12 +
//                        static_cast<ScalarType>(28) * Y22)) +
//            THREE * (-TWO * BB *
//                         (static_cast<ScalarType>(46) * Y02 * Y11 -
//                          static_cast<ScalarType>(14) * Y0 * Y111 + Y1111 -
//                          static_cast<ScalarType>(14) * Y111 * Y2 -
//                          static_cast<ScalarType>(14) * Y01 * Y22 +
//                          static_cast<ScalarType>(28) * Y11 * Y22 +
//                          Y00 * (static_cast<ScalarType>(28) * Y11 -
//                                 static_cast<ScalarType>(14) * Y12 + Y22)) +
//                     FIVE * Z00 + static_cast<ScalarType>(40) * Z01p -
//                     TWO * Z02 + static_cast<ScalarType>(8) * Z1p1p +
//                     static_cast<ScalarType>(40) * Z1p2 + FIVE * Z22),
//        -TWO * AA *
//                (static_cast<ScalarType>(3) * X02 * X11 -
//                 static_cast<ScalarType>(7) * X0 * X111 + THREE * X1111 -
//                 static_cast<ScalarType>(7) * X111 * X2) -
//            SIX * BB * Y02 * Y11 + static_cast<ScalarType>(14) * BB * Y0 *
//            Y111 - SIX * BB * Y1111 - TWO * AB *
//                (static_cast<ScalarType>(2) * X12 * Y01 -
//                 static_cast<ScalarType>(7) * X01 * Y11 + X02 * Y11 -
//                 static_cast<ScalarType>(7) * X12 * Y11 +
//                 X11 * (-static_cast<ScalarType>(7) * Y01 + Y02 + SIX * Y11
//                 -
//                        static_cast<ScalarType>(7) * Y12) +
//                 TWO * X01 * Y12) +
//            static_cast<ScalarType>(14) * BB * Y111 * Y2 - FIVE * Z01p + Z02
//            - static_cast<ScalarType>(7) * Z1p1p - FIVE * Z1p2,
//        -(AA * X1111) - TWO * AB * X11 * Y11 - BB * Y1111 + Z1p1p});
//   for (size_t i = 0; i < 3; ++i) {
//     moments.volume() += coeffs[i] * m0_basis[i];
//   }
//   for (size_t i = 0; i < 4; ++i) {
//     moments.centroid()[0] += coeffs[3 + i] * m1x_basis[i];
//     moments.centroid()[1] += coeffs[3 + i] * m1y_basis[i];
//   }
//   for (size_t i = 0; i < 5; ++i) {
//     moments.centroid()[2] += coeffs[7 + i] * m1z_basis[i];
//   }
//   moments.volume() *= area_proj_triangle;
//   moments.centroid()[0] *= area_proj_triangle;
//   moments.centroid()[1] *= area_proj_triangle;
//   moments.centroid()[2] *= area_proj_triangle;
//   return moments;
// }  // namespace IRL

template <class ReturnType, class ScalarType>
ReturnType computeFaceOnlyContribution(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PlaneBase<ScalarType>& a_face_plane,
    const PtBase<ScalarType>& a_pt_ref) {
  if constexpr (std::is_same_v<ReturnType, VolumeBase<double>> ||
                std::is_same_v<ReturnType, VolumeBase<Quad_t>>) {
    /* Defining constants and types */
    const ScalarType EPSILON = machine_epsilon<ScalarType>();
    const ScalarType ZERO = static_cast<ScalarType>(0);
    const ScalarType FOUR = static_cast<ScalarType>(4);

    /* Function */
    assert(a_paraboloid.a() * a_paraboloid.b() > ZERO);
    assert(fabs(a_face_plane.normal()[2]) > EPSILON);
    const ScalarType a =
        -a_face_plane.normal()[0] / safelyEpsilon(a_face_plane.normal()[2]);
    const ScalarType b =
        -a_face_plane.normal()[1] / safelyEpsilon(a_face_plane.normal()[2]);
    const ScalarType c =
        a_face_plane.distance() / safelyEpsilon(a_face_plane.normal()[2]);
    const ScalarType factor = FOUR * a_paraboloid.a() * a_paraboloid.b() * c -
                              a_paraboloid.a() * b * b -
                              a_paraboloid.b() * a * a;
    return copysign(machine_pi<ScalarType>() * factor * factor /
                        (static_cast<ScalarType>(32) *
                         pow(a_paraboloid.a() * a_paraboloid.b(),
                             static_cast<ScalarType>(2.5))),
                    -a_face_plane.normal()[2]);
  } else if constexpr (std::is_same_v<ReturnType, VolumeMomentsBase<double>> ||
                       std::is_same_v<ReturnType, VolumeMomentsBase<Quad_t>>) {
    /* Defining constants and types */
    const ScalarType EPSILON = machine_epsilon<ScalarType>();
    const ScalarType ZERO = static_cast<ScalarType>(0);
    const ScalarType FOUR = static_cast<ScalarType>(4);
    const ScalarType FIVE = static_cast<ScalarType>(5);

    /* Function */
    assert(a_paraboloid.a() * a_paraboloid.b() > ZERO);
    assert(fabs(a_face_plane.normal()[2]) > EPSILON);
    const ScalarType a = -a_face_plane.normal()[0] / a_face_plane.normal()[2];
    const ScalarType b = -a_face_plane.normal()[1] / a_face_plane.normal()[2];
    const ScalarType c = a_face_plane.distance() / a_face_plane.normal()[2];
    const auto A = a_paraboloid.a(), B = a_paraboloid.b();
    auto moments = ReturnType::fromScalarConstant(ZERO);
    const ScalarType factor = (a * a * B + A * (b * b - FOUR * B * c)) *
                              (a * a * B + A * (b * b - FOUR * B * c)) *
                              machine_pi<ScalarType>();
    moments.volume() =
        copysign(factor / (static_cast<ScalarType>(32) *
                           pow(A * B, static_cast<ScalarType>(2.5))),
                 -a_face_plane.normal()[2]);
    moments.centroid()[0] =
        a * B *
        copysign(factor / (static_cast<ScalarType>(64) *
                           pow(A * B, static_cast<ScalarType>(3.5))),
                 a_face_plane.normal()[2]);
    moments.centroid()[1] =
        b * A *
        copysign(factor / (static_cast<ScalarType>(64) *
                           pow(A * B, static_cast<ScalarType>(3.5))),
                 a_face_plane.normal()[2]);
    moments.centroid()[2] =
        (FIVE * A * (b * b) + FIVE * (a * a) * B -
         static_cast<ScalarType>(8) * A * B * c) *
        copysign(factor / (static_cast<ScalarType>(384) *
                           pow(A * B, static_cast<ScalarType>(3.5))),
                 a_face_plane.normal()[2]);
    return moments;
  }
}

template <class ReturnType, class ScalarType>
inline ReturnType computeFaceOnlyContributionWithGradient(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PlaneBase<ScalarType>& a_face_plane,
    const PtWithGradientBase<typename ReturnType::gradient_type, ScalarType>&
        a_pt_ref) {
  /* Defining constants and types */
  using gradient_type = typename ReturnType::gradient_type;
  using scalar_container_type = ScalarWithGradient<gradient_type>;
  const ScalarType EPSILON = machine_epsilon<ScalarType>();
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType FOUR = static_cast<ScalarType>(4);
  const ScalarType FIVE = static_cast<ScalarType>(5);

  /* Function */
  auto moments_with_gradient = ReturnType::fromScalarConstant(ZERO);
  assert(a_paraboloid.a() * a_paraboloid.b() > ZERO);
  assert(fabs(a_face_plane.normal()[2]) > EPSILON);

  // Compute plane gradient
  const auto& pt_ref = a_pt_ref.getPt();
  const auto& pt_ref_grad = a_pt_ref.getData();
  auto face_normal = a_face_plane.normal();
  auto face_normal_withgrad = PtWithGradient<gradient_type>(
      PtBase<ScalarType>(face_normal[0], face_normal[1], face_normal[2]));
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
  auto face_distance = face_normal[0] * pt_ref[0] + face_normal[1] * pt_ref[1] +
                       face_normal[2] * pt_ref[2];
  auto face_distance_grad =
      face_normal_grad[0] * pt_ref[0] + face_normal_grad[1] * pt_ref[1] +
      face_normal_grad[2] * pt_ref[2] + face_normal[0] * pt_ref_grad[0] +
      face_normal[1] * pt_ref_grad[1] + face_normal[2] * pt_ref_grad[2];

  const ScalarType a_value = -face_normal[0] / safelyEpsilon(face_normal[2]);
  const ScalarType b_value = -face_normal[1] / safelyEpsilon(face_normal[2]);
  const ScalarType c_value = face_distance / safelyEpsilon(face_normal[2]);
  const auto a_grad = -face_normal_grad[0] / safelyEpsilon(face_normal[2]) +
                      face_normal[0] * face_normal_grad[2] /
                          safelyEpsilon(face_normal[2] * face_normal[2]);
  const auto b_grad = -face_normal_grad[1] / safelyEpsilon(face_normal[2]) +
                      face_normal[1] * face_normal_grad[2] /
                          safelyEpsilon(face_normal[2] * face_normal[2]);
  const auto c_grad = face_distance_grad / safelyEpsilon(face_normal[2]) -
                      face_distance * face_normal_grad[2] /
                          safelyEpsilon(face_normal[2] * face_normal[2]);
  const auto a = scalar_container_type(a_value, a_grad);
  const auto b = scalar_container_type(b_value, b_grad);
  const auto c = scalar_container_type(c_value, c_grad);
  auto A_grad = gradient_type(ZERO), B_grad = gradient_type(ZERO);
  A_grad.setGradA(ONE);
  B_grad.setGradB(ONE);
  const auto A = scalar_container_type(a_paraboloid.a(), A_grad),
             B = scalar_container_type(a_paraboloid.b(), B_grad);
  const auto factor = (a * a * B + A * (b * b - FOUR * B * c)) *
                      (a * a * B + A * (b * b - FOUR * B * c)) *
                      machine_pi<ScalarType>();
  const ScalarType sign_normal = copysign(ONE, a_face_plane.normal()[2]);
  const auto m0 =
      -sign_normal * factor /
      (static_cast<ScalarType>(32) *
       PowMoments<scalar_container_type>(A * B, static_cast<ScalarType>(2.5)));
  moments_with_gradient.volume() = m0.value();
  moments_with_gradient.volume_gradient() = m0.gradient();
  if constexpr (std::is_same<typename ReturnType::moment_type,
                             VolumeMoments>::value) {
    const auto m1x = sign_normal * a * B * factor /
                     (static_cast<ScalarType>(64) *
                      PowMoments<scalar_container_type>(
                          A * B, static_cast<ScalarType>(3.5)));
    const auto m1y = sign_normal * b * A * factor /
                     (static_cast<ScalarType>(64) *
                      PowMoments<scalar_container_type>(
                          A * B, static_cast<ScalarType>(3.5)));
    const auto m1z = sign_normal *
                     (FIVE * A * (b * b) + FIVE * (a * a) * B -
                      static_cast<ScalarType>(8) * A * B * c) *
                     factor /
                     (static_cast<ScalarType>(384) *
                      PowMoments<scalar_container_type>(
                          A * B, static_cast<ScalarType>(3.5)));
    PtBase<ScalarType>& centroid = moments_with_gradient.centroid().getPt();
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

template <class ReturnType, class ScalarType>
ReturnType computeTriangleCorrection(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1,
    const PtBase<ScalarType>& a_pt_2) {
  if constexpr (std::is_same_v<ReturnType, VolumeBase<double>> ||
                std::is_same_v<ReturnType, VolumeBase<Quad_t>>) {
    return (-a_paraboloid.a() * (a_pt_0[0] + a_pt_1[0]) *
                (a_pt_1[0] + a_pt_2[0]) +
            -a_paraboloid.b() * (a_pt_0[1] + a_pt_1[1]) *
                (a_pt_1[1] + a_pt_2[1]) -
            a_pt_0[2] - static_cast<ScalarType>(2) * a_pt_1[2] - a_pt_2[2]) /
           static_cast<ScalarType>(12) *
           ((a_pt_1[1] - a_pt_2[1]) * a_pt_0[0] +
            (a_pt_2[1] - a_pt_0[1]) * a_pt_1[0] +
            (a_pt_0[1] - a_pt_1[1]) * a_pt_2[0]);
  } else if constexpr (std::is_same_v<ReturnType, VolumeMomentsBase<double>> ||
                       std::is_same_v<ReturnType, VolumeMomentsBase<Quad_t>>) {
    /* Defining constants and types */
    const ScalarType ZERO = static_cast<ScalarType>(0);
    const ScalarType ONE = static_cast<ScalarType>(1);
    const ScalarType TWO = static_cast<ScalarType>(2);
    const ScalarType THREE = static_cast<ScalarType>(3);
    const ScalarType FOUR = static_cast<ScalarType>(4);
    const ScalarType FIVE = static_cast<ScalarType>(5);
    const ScalarType SIX = static_cast<ScalarType>(6);
    const ScalarType HALF = ONE / TWO;

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ZERO);
    const ScalarType A = a_paraboloid.a(), B = a_paraboloid.b();
    const ScalarType X0 = a_pt_0[0], X1 = a_pt_1[0], X2 = a_pt_2[0];
    const ScalarType Y0 = a_pt_0[1], Y1 = a_pt_1[1], Y2 = a_pt_2[1];
    const ScalarType Z0 = a_pt_0[2], Z1 = a_pt_1[2], Z2 = a_pt_2[2];
    const ScalarType triangle_area =
        HALF * ((a_pt_1[1] - a_pt_2[1]) * a_pt_0[0] +
                (a_pt_2[1] - a_pt_0[1]) * a_pt_1[0] +
                (a_pt_0[1] - a_pt_1[1]) * a_pt_2[0]);
    moments.volume() =
        (-a_paraboloid.a() * (a_pt_0[0] + a_pt_1[0]) * (a_pt_1[0] + a_pt_2[0]) +
         -a_paraboloid.b() * (a_pt_0[1] + a_pt_1[1]) * (a_pt_1[1] + a_pt_2[1]) -
         a_pt_0[2] - TWO * a_pt_1[2] - a_pt_2[2]) *
        triangle_area / SIX;
    moments.centroid()[0] =
        triangle_area *
        ((A * (-(X0 * X0 * X0) - X1 * X1 * X1 - X1 * X1 * X2 - X1 * (X2 * X2) -
               X2 * X2 * X2 - X0 * X0 * (X1 + X2) -
               X0 * (X1 * X1 + X1 * X2 + X2 * X2))) /
             static_cast<ScalarType>(10) +
         (B * (-(X1 * (Y0 * Y0 + TWO * Y0 * Y1 + THREE * (Y1 * Y1) + Y0 * Y2 +
                       TWO * Y1 * Y2 + Y2 * Y2)) -
               X2 * (Y0 * Y0 + Y0 * Y1 + Y1 * Y1 + TWO * Y0 * Y2 +
                     TWO * Y1 * Y2 + THREE * (Y2 * Y2)) -
               X0 * (THREE * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
                     TWO * Y0 * (Y1 + Y2)))) /
             static_cast<ScalarType>(30) +
         (-(X0 * (TWO * Z0 + Z1 + Z2)) - X1 * (Z0 + TWO * Z1 + Z2) -
          X2 * (Z0 + Z1 + TWO * Z2)) /
             static_cast<ScalarType>(12));
    moments.centroid()[1] =
        -triangle_area *
        ((B * (Y0 * Y0 * Y0 + Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) +
               Y2 * Y2 * Y2 + Y0 * Y0 * (Y1 + Y2) +
               Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2))) /
             static_cast<ScalarType>(10) +
         (A *
          (X0 * X0 * (THREE * Y0 + Y1 + Y2) + X1 * X1 * (Y0 + THREE * Y1 + Y2) +
           X2 * X2 * (Y0 + Y1 + THREE * Y2) + X1 * X2 * (Y0 + TWO * (Y1 + Y2)) +
           X0 * (X1 * (TWO * Y0 + TWO * Y1 + Y2) +
                 X2 * (TWO * Y0 + Y1 + TWO * Y2)))) /
             static_cast<ScalarType>(30) +
         (Y0 * (TWO * Z0 + Z1 + Z2) + Y1 * (Z0 + TWO * Z1 + Z2) +
          Y2 * (Z0 + Z1 + TWO * Z2)) /
             static_cast<ScalarType>(12));
    moments.centroid()[2] =
        triangle_area *
        ((A * A *
          (X0 * X0 * X0 * X0 + X1 * X1 * X1 * X1 + X1 * X1 * X1 * X2 +
           X1 * X1 * (X2 * X2) + X1 * (X2 * X2 * X2) + X2 * X2 * X2 * X2 +
           X0 * X0 * X0 * (X1 + X2) + X0 * X0 * (X1 * X1 + X1 * X2 + X2 * X2) +
           X0 *
               (X1 * X1 * X1 + X1 * X1 * X2 + X1 * (X2 * X2) + X2 * X2 * X2))) /
             static_cast<ScalarType>(30) +
         (B * B *
          (Y0 * Y0 * Y0 * Y0 + Y1 * Y1 * Y1 * Y1 + Y1 * Y1 * Y1 * Y2 +
           Y1 * Y1 * (Y2 * Y2) + Y1 * (Y2 * Y2 * Y2) + Y2 * Y2 * Y2 * Y2 +
           Y0 * Y0 * Y0 * (Y1 + Y2) + Y0 * Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2) +
           Y0 *
               (Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) + Y2 * Y2 * Y2))) /
             static_cast<ScalarType>(30) +
         (A * B *
          (X1 * X2 *
               (Y0 * Y0 + THREE * (Y1 * Y1) + FOUR * Y1 * Y2 +
                THREE * (Y2 * Y2) + TWO * Y0 * (Y1 + Y2)) +
           X0 * X0 *
               (SIX * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
                THREE * Y0 * (Y1 + Y2)) +
           X1 * X1 *
               (Y0 * Y0 + SIX * (Y1 * Y1) + THREE * Y1 * Y2 + Y2 * Y2 +
                Y0 * (THREE * Y1 + Y2)) +
           X2 * X2 *
               (Y0 * Y0 + Y1 * Y1 + THREE * Y1 * Y2 + SIX * (Y2 * Y2) +
                Y0 * (Y1 + THREE * Y2)) +
           X0 * (X1 * (THREE * (Y0 * Y0) + FOUR * Y0 * Y1 + THREE * (Y1 * Y1) +
                       TWO * Y0 * Y2 + TWO * Y1 * Y2 + Y2 * Y2) +
                 X2 * (THREE * (Y0 * Y0) + TWO * Y0 * Y1 + Y1 * Y1 +
                       FOUR * Y0 * Y2 + TWO * Y1 * Y2 + THREE * (Y2 * Y2))))) /
             static_cast<ScalarType>(90) +
         (-(Z0 * Z0) - Z1 * Z1 - Z1 * Z2 - Z2 * Z2 - Z0 * (Z1 + Z2)) /
             static_cast<ScalarType>(12));
    return moments;
  }
}

template <class ReturnType, class ScalarType>
inline ReturnType computeTriangleCorrectionWithGradient(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtWithGradientBase<typename ReturnType::gradient_type, ScalarType>&
        a_pt_0,
    const PtWithGradientBase<typename ReturnType::gradient_type, ScalarType>&
        a_pt_1,
    const PtWithGradientBase<typename ReturnType::gradient_type, ScalarType>&
        a_pt_2) {
  /* Defining constants and types */
  using gradient_type = typename ReturnType::gradient_type;
  using scalar_container_type = ScalarWithGradient<gradient_type>;
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType THREE = static_cast<ScalarType>(3);
  const ScalarType FOUR = static_cast<ScalarType>(4);
  const ScalarType FIVE = static_cast<ScalarType>(5);
  const ScalarType SIX = static_cast<ScalarType>(6);
  const ScalarType HALF = ONE / TWO;

  /* Function */
  auto moments_with_gradient = ReturnType::fromScalarConstant(ZERO);
  std::array<scalar_container_type, 3> pt_0, pt_1, pt_2;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    pt_0[d] = scalar_container_type(a_pt_0.getPt()[d], a_pt_0.getData()[d]);
    pt_1[d] = scalar_container_type(a_pt_1.getPt()[d], a_pt_1.getData()[d]);
    pt_2[d] = scalar_container_type(a_pt_2.getPt()[d], a_pt_2.getData()[d]);
  }
  auto A_grad = gradient_type(ZERO), B_grad = gradient_type(ZERO);
  A_grad.setGradA(ONE);
  B_grad.setGradB(ONE);
  const auto A = scalar_container_type(a_paraboloid.a(), A_grad),
             B = scalar_container_type(a_paraboloid.b(), B_grad);
  const auto X0 = pt_0[0], X1 = pt_1[0], X2 = pt_2[0];
  const auto Y0 = pt_0[1], Y1 = pt_1[1], Y2 = pt_2[1];
  const auto Z0 = pt_0[2], Z1 = pt_1[2], Z2 = pt_2[2];
  const auto triangle_area =
      HALF * ((pt_1[1] - pt_2[1]) * pt_0[0] + (pt_2[1] - pt_0[1]) * pt_1[0] +
              (pt_0[1] - pt_1[1]) * pt_2[0]);
  const auto m0 = (-A * (pt_0[0] + pt_1[0]) * (pt_1[0] + pt_2[0]) +
                   -B * (pt_0[1] + pt_1[1]) * (pt_1[1] + pt_2[1]) - pt_0[2] -
                   TWO * pt_1[2] - pt_2[2]) *
                  triangle_area / SIX;
  moments_with_gradient.volume() = m0.value();
  moments_with_gradient.volume_gradient() = m0.gradient();
  if constexpr (std::is_same<typename ReturnType::moment_type,
                             VolumeMoments>::value) {
    const auto m1x =
        triangle_area *
        ((A * (-(X0 * X0 * X0) - X1 * X1 * X1 - X1 * X1 * X2 - X1 * (X2 * X2) -
               X2 * X2 * X2 - X0 * X0 * (X1 + X2) -
               X0 * (X1 * X1 + X1 * X2 + X2 * X2))) /
             static_cast<ScalarType>(10) +
         (B * (-(X1 * (Y0 * Y0 + TWO * Y0 * Y1 + THREE * (Y1 * Y1) + Y0 * Y2 +
                       TWO * Y1 * Y2 + Y2 * Y2)) -
               X2 * (Y0 * Y0 + Y0 * Y1 + Y1 * Y1 + TWO * Y0 * Y2 +
                     TWO * Y1 * Y2 + THREE * (Y2 * Y2)) -
               X0 * (THREE * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
                     TWO * Y0 * (Y1 + Y2)))) /
             static_cast<ScalarType>(30) +
         (-(X0 * (TWO * Z0 + Z1 + Z2)) - X1 * (Z0 + TWO * Z1 + Z2) -
          X2 * (Z0 + Z1 + TWO * Z2)) /
             static_cast<ScalarType>(12));
    const auto m1y =
        -triangle_area *
        ((B * (Y0 * Y0 * Y0 + Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) +
               Y2 * Y2 * Y2 + Y0 * Y0 * (Y1 + Y2) +
               Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2))) /
             static_cast<ScalarType>(10) +
         (A *
          (X0 * X0 * (THREE * Y0 + Y1 + Y2) + X1 * X1 * (Y0 + THREE * Y1 + Y2) +
           X2 * X2 * (Y0 + Y1 + THREE * Y2) + X1 * X2 * (Y0 + TWO * (Y1 + Y2)) +
           X0 * (X1 * (TWO * Y0 + TWO * Y1 + Y2) +
                 X2 * (TWO * Y0 + Y1 + TWO * Y2)))) /
             static_cast<ScalarType>(30) +
         (Y0 * (TWO * Z0 + Z1 + Z2) + Y1 * (Z0 + TWO * Z1 + Z2) +
          Y2 * (Z0 + Z1 + TWO * Z2)) /
             static_cast<ScalarType>(12));
    const auto m1z =
        triangle_area *
        ((A * A *
          (X0 * X0 * X0 * X0 + X1 * X1 * X1 * X1 + X1 * X1 * X1 * X2 +
           X1 * X1 * (X2 * X2) + X1 * (X2 * X2 * X2) + X2 * X2 * X2 * X2 +
           X0 * X0 * X0 * (X1 + X2) + X0 * X0 * (X1 * X1 + X1 * X2 + X2 * X2) +
           X0 *
               (X1 * X1 * X1 + X1 * X1 * X2 + X1 * (X2 * X2) + X2 * X2 * X2))) /
             static_cast<ScalarType>(30) +
         (B * B *
          (Y0 * Y0 * Y0 * Y0 + Y1 * Y1 * Y1 * Y1 + Y1 * Y1 * Y1 * Y2 +
           Y1 * Y1 * (Y2 * Y2) + Y1 * (Y2 * Y2 * Y2) + Y2 * Y2 * Y2 * Y2 +
           Y0 * Y0 * Y0 * (Y1 + Y2) + Y0 * Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2) +
           Y0 *
               (Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) + Y2 * Y2 * Y2))) /
             static_cast<ScalarType>(30) +
         (A * B *
          (X1 * X2 *
               (Y0 * Y0 + THREE * (Y1 * Y1) + FOUR * Y1 * Y2 +
                THREE * (Y2 * Y2) + TWO * Y0 * (Y1 + Y2)) +
           X0 * X0 *
               (SIX * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
                THREE * Y0 * (Y1 + Y2)) +
           X1 * X1 *
               (Y0 * Y0 + SIX * (Y1 * Y1) + THREE * Y1 * Y2 + Y2 * Y2 +
                Y0 * (THREE * Y1 + Y2)) +
           X2 * X2 *
               (Y0 * Y0 + Y1 * Y1 + THREE * Y1 * Y2 + SIX * (Y2 * Y2) +
                Y0 * (Y1 + THREE * Y2)) +
           X0 * (X1 * (THREE * (Y0 * Y0) + FOUR * Y0 * Y1 + THREE * (Y1 * Y1) +
                       TWO * Y0 * Y2 + TWO * Y1 * Y2 + Y2 * Y2) +
                 X2 * (THREE * (Y0 * Y0) + TWO * Y0 * Y1 + Y1 * Y1 +
                       FOUR * Y0 * Y2 + TWO * Y1 * Y2 + THREE * (Y2 * Y2))))) /
             static_cast<ScalarType>(90) +
         (-(Z0 * Z0) - Z1 * Z1 - Z1 * Z2 - Z2 * Z2 - Z0 * (Z1 + Z2)) /
             static_cast<ScalarType>(12));
    PtBase<ScalarType>& centroid = moments_with_gradient.centroid().getPt();
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

}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_MOMENT_CONTRIBUTIONS_TPP_
