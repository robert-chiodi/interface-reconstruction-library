// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_AMR_TPP_
#define IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_AMR_TPP_

#include <float.h>
#include <cassert>
#include <cmath>
#include <fstream>

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"

namespace IRL {

const double triangleSignedArea(const Pt& a_pt_0, const Pt& a_pt_1,
                                const Pt& a_pt_2) {
  return 0.5 * ((a_pt_1[0] - a_pt_0[0]) * (a_pt_2[1] - a_pt_0[1]) -
                (a_pt_2[0] - a_pt_0[0]) * (a_pt_1[1] - a_pt_0[1]));
}

template <>
void kahanSummationMoments(
    std::array<std::pair<Volume, Volume>, N_AMR_STRATEGIES>& a_full_moments,
    std::array<std::pair<Volume, Volume>, N_AMR_STRATEGIES>& a_full_moments_ref,
    std::array<Volume, N_AMR_STRATEGIES>& a_moments_to_add) {
  for (UnsignedIndex_t m = 0; m < N_AMR_STRATEGIES; ++m) {
    const double y_first = a_moments_to_add[m] - a_full_moments_ref[m].first;
    const double y_second =
        std::fabs(a_moments_to_add[m]) - a_full_moments_ref[m].second;
    const double t_first = a_full_moments[m].first + y_first;
    const double t_second = a_full_moments[m].first + y_second;
    a_full_moments_ref[m].first = (t_first - a_full_moments[m].first) - y_first;
    a_full_moments_ref[m].second =
        (t_second - a_full_moments[m].second) - y_second;
    a_full_moments[m].first = t_first;
    a_full_moments[m].second = t_second;
  }
}

template <>
void kahanSummationMoments(
    std::array<std::pair<VolumeMoments, VolumeMoments>, N_AMR_STRATEGIES>&
        a_full_moments,
    std::array<std::pair<VolumeMoments, VolumeMoments>, N_AMR_STRATEGIES>&
        a_full_moments_ref,
    std::array<VolumeMoments, N_AMR_STRATEGIES>& a_moments_to_add) {
  for (UnsignedIndex_t m = 0; m < N_AMR_STRATEGIES; ++m) {
    const double y_first_m0 =
        a_moments_to_add[m].volume() - a_full_moments_ref[m].first.volume();
    const double y_first_m1x = a_moments_to_add[m].centroid()[0] -
                               a_full_moments_ref[m].first.centroid()[0];
    const double y_first_m1y = a_moments_to_add[m].centroid()[1] -
                               a_full_moments_ref[m].first.centroid()[1];
    const double y_first_m1z = a_moments_to_add[m].centroid()[2] -
                               a_full_moments_ref[m].first.centroid()[2];
    const double y_second_m0 = std::fabs(a_moments_to_add[m].volume()) -
                               a_full_moments_ref[m].second.volume();
    const double y_second_m1x = std::fabs(a_moments_to_add[m].centroid()[0]) -
                                a_full_moments_ref[m].second.centroid()[0];
    const double y_second_m1y = std::fabs(a_moments_to_add[m].centroid()[1]) -
                                a_full_moments_ref[m].second.centroid()[1];
    const double y_second_m1z = std::fabs(a_moments_to_add[m].centroid()[2]) -
                                a_full_moments_ref[m].second.centroid()[2];
    const double t_first_m0 = a_full_moments[m].first.volume() + y_first_m0;
    const double t_first_m1x =
        a_full_moments[m].first.centroid()[0] + y_first_m1x;
    const double t_first_m1y =
        a_full_moments[m].first.centroid()[1] + y_first_m1y;
    const double t_first_m1z =
        a_full_moments[m].first.centroid()[2] + y_first_m1z;
    const double t_second_m0 = a_full_moments[m].first.volume() + y_second_m0;
    const double t_second_m1x =
        a_full_moments[m].first.centroid()[0] + y_second_m1x;
    const double t_second_m1y =
        a_full_moments[m].first.centroid()[1] + y_second_m1y;
    const double t_second_m1z =
        a_full_moments[m].first.centroid()[2] + y_second_m1z;
    a_full_moments_ref[m].first.volume() =
        (t_first_m0 - a_full_moments[m].first.volume()) - y_first_m0;
    a_full_moments_ref[m].first.centroid()[0] =
        (t_first_m1x - a_full_moments[m].first.centroid()[0]) - y_first_m1x;
    a_full_moments_ref[m].first.centroid()[1] =
        (t_first_m1y - a_full_moments[m].first.centroid()[1]) - y_first_m1y;
    a_full_moments_ref[m].first.centroid()[2] =
        (t_first_m1z - a_full_moments[m].first.centroid()[2]) - y_first_m1z;
    a_full_moments_ref[m].second.volume() =
        (t_second_m0 - a_full_moments[m].second.volume()) - y_second_m0;
    a_full_moments_ref[m].second.centroid()[0] =
        (t_second_m1x - a_full_moments[m].second.centroid()[0]) - y_second_m1x;
    a_full_moments_ref[m].second.centroid()[1] =
        (t_second_m1y - a_full_moments[m].second.centroid()[1]) - y_second_m1y;
    a_full_moments_ref[m].second.centroid()[2] =
        (t_second_m1z - a_full_moments[m].second.centroid()[2]) - y_second_m1z;
    a_full_moments[m].first.volume() = t_first_m0;
    a_full_moments[m].first.centroid()[0] = t_first_m1x;
    a_full_moments[m].first.centroid()[1] = t_first_m1y;
    a_full_moments[m].first.centroid()[2] = t_first_m1z;
    a_full_moments[m].second.volume() = t_second_m0;
    a_full_moments[m].second.centroid()[0] = t_second_m1x;
    a_full_moments[m].second.centroid()[1] = t_second_m1y;
    a_full_moments[m].second.centroid()[2] = t_second_m1z;
  }
}

template <>
inline Volume computeMomentContributionClippedTriangle<Volume, 0>(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print) {
  if (a_print) {
    amr_triangles_clipped.insert(
        amr_triangles_clipped.end(),
        {a_pt_0[0], a_pt_0[1], a_pt_0[2], a_pt_1[0], a_pt_1[1], a_pt_1[2],
         a_pt_2[0], a_pt_2[1], a_pt_2[2]});
  }
  return -(a_aligned_paraboloid.a() *
               (a_pt_0[0] * (a_pt_0[0] + a_pt_1[0] + a_pt_2[0]) +
                a_pt_1[0] * (a_pt_1[0] + a_pt_2[0]) + a_pt_2[0] * a_pt_2[0]) +
           a_aligned_paraboloid.b() *
               (a_pt_0[1] * (a_pt_0[1] + a_pt_1[1] + a_pt_2[1]) +
                a_pt_1[1] * (a_pt_1[1] + a_pt_2[1]) + a_pt_2[1] * a_pt_2[1])) *
         a_signed_area / 6.0;
}

template <>
inline VolumeMoments computeMomentContributionClippedTriangle<VolumeMoments, 0>(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print) {
  if (a_print) {
    amr_triangles_clipped.insert(
        amr_triangles_clipped.end(),
        {a_pt_0[0], a_pt_0[1], a_pt_0[2], a_pt_1[0], a_pt_1[1], a_pt_1[2],
         a_pt_2[0], a_pt_2[1], a_pt_2[2]});
  }
  auto moments = VolumeMoments::fromScalarConstant(0.0);
  moments.volume() =
      -(a_aligned_paraboloid.a() *
            (a_pt_0[0] * (a_pt_0[0] + a_pt_1[0] + a_pt_2[0]) +
             a_pt_1[0] * (a_pt_1[0] + a_pt_2[0]) + a_pt_2[0] * a_pt_2[0]) +
        a_aligned_paraboloid.b() *
            (a_pt_0[1] * (a_pt_0[1] + a_pt_1[1] + a_pt_2[1]) +
             a_pt_1[1] * (a_pt_1[1] + a_pt_2[1]) + a_pt_2[1] * a_pt_2[1])) *
      a_signed_area / 6.0;
  moments.centroid()[0] =
      -(a_aligned_paraboloid.a() * (6. * (a_pt_0[0] * a_pt_0[0] * a_pt_0[0]) +
                                    6. * (a_pt_0[0] * a_pt_0[0]) * a_pt_1[0] +
                                    6. * a_pt_0[0] * (a_pt_1[0] * a_pt_1[0]) +
                                    6. * (a_pt_1[0] * a_pt_1[0] * a_pt_1[0]) +
                                    6. * (a_pt_0[0] * a_pt_0[0]) * a_pt_2[0] +
                                    6. * a_pt_0[0] * a_pt_1[0] * a_pt_2[0] +
                                    6. * (a_pt_1[0] * a_pt_1[0]) * a_pt_2[0] +
                                    6. * a_pt_0[0] * (a_pt_2[0] * a_pt_2[0]) +
                                    6. * a_pt_1[0] * (a_pt_2[0] * a_pt_2[0]) +
                                    6. * (a_pt_2[0] * a_pt_2[0] * a_pt_2[0])) +
        a_aligned_paraboloid.b() * (6. * a_pt_0[0] * (a_pt_0[1] * a_pt_0[1]) +
                                    2. * a_pt_1[0] * (a_pt_0[1] * a_pt_0[1]) +
                                    2. * a_pt_2[0] * (a_pt_0[1] * a_pt_0[1]) +
                                    4. * a_pt_0[0] * a_pt_0[1] * a_pt_1[1] +
                                    4. * a_pt_1[0] * a_pt_0[1] * a_pt_1[1] +
                                    2. * a_pt_2[0] * a_pt_0[1] * a_pt_1[1] +
                                    2. * a_pt_0[0] * (a_pt_1[1] * a_pt_1[1]) +
                                    6. * a_pt_1[0] * (a_pt_1[1] * a_pt_1[1]) +
                                    2. * a_pt_2[0] * (a_pt_1[1] * a_pt_1[1]) +
                                    4. * a_pt_0[0] * a_pt_0[1] * a_pt_2[1] +
                                    2. * a_pt_1[0] * a_pt_0[1] * a_pt_2[1] +
                                    4. * a_pt_2[0] * a_pt_0[1] * a_pt_2[1] +
                                    2. * a_pt_0[0] * a_pt_1[1] * a_pt_2[1] +
                                    4. * a_pt_1[0] * a_pt_1[1] * a_pt_2[1] +
                                    4. * a_pt_2[0] * a_pt_1[1] * a_pt_2[1] +
                                    2. * a_pt_0[0] * (a_pt_2[1] * a_pt_2[1]) +
                                    2. * a_pt_1[0] * (a_pt_2[1] * a_pt_2[1]) +
                                    6. * a_pt_2[0] * (a_pt_2[1] * a_pt_2[1]))) *
      a_signed_area / 60.0;
  moments.centroid()[1] =
      -(a_aligned_paraboloid.b() * (6. * (a_pt_0[1] * a_pt_0[1] * a_pt_0[1]) +
                                    6. * (a_pt_0[1] * a_pt_0[1]) * a_pt_1[1] +
                                    6. * a_pt_0[1] * (a_pt_1[1] * a_pt_1[1]) +
                                    6. * (a_pt_1[1] * a_pt_1[1] * a_pt_1[1]) +
                                    6. * (a_pt_0[1] * a_pt_0[1]) * a_pt_2[1] +
                                    6. * a_pt_0[1] * a_pt_1[1] * a_pt_2[1] +
                                    6. * (a_pt_1[1] * a_pt_1[1]) * a_pt_2[1] +
                                    6. * a_pt_0[1] * (a_pt_2[1] * a_pt_2[1]) +
                                    6. * a_pt_1[1] * (a_pt_2[1] * a_pt_2[1]) +
                                    6. * (a_pt_2[1] * a_pt_2[1] * a_pt_2[1])) +
        a_aligned_paraboloid.a() * (6. * a_pt_0[1] * (a_pt_0[0] * a_pt_0[0]) +
                                    2. * a_pt_1[1] * (a_pt_0[0] * a_pt_0[0]) +
                                    2. * a_pt_2[1] * (a_pt_0[0] * a_pt_0[0]) +
                                    4. * a_pt_0[1] * a_pt_0[0] * a_pt_1[0] +
                                    4. * a_pt_1[1] * a_pt_0[0] * a_pt_1[0] +
                                    2. * a_pt_2[1] * a_pt_0[0] * a_pt_1[0] +
                                    2. * a_pt_0[1] * (a_pt_1[0] * a_pt_1[0]) +
                                    6. * a_pt_1[1] * (a_pt_1[0] * a_pt_1[0]) +
                                    2. * a_pt_2[1] * (a_pt_1[0] * a_pt_1[0]) +
                                    4. * a_pt_0[1] * a_pt_0[0] * a_pt_2[0] +
                                    2. * a_pt_1[1] * a_pt_0[0] * a_pt_2[0] +
                                    4. * a_pt_2[1] * a_pt_0[0] * a_pt_2[0] +
                                    2. * a_pt_0[1] * a_pt_1[0] * a_pt_2[0] +
                                    4. * a_pt_1[1] * a_pt_1[0] * a_pt_2[0] +
                                    4. * a_pt_2[1] * a_pt_1[0] * a_pt_2[0] +
                                    2. * a_pt_0[1] * (a_pt_2[0] * a_pt_2[0]) +
                                    2. * a_pt_1[1] * (a_pt_2[0] * a_pt_2[0]) +
                                    6. * a_pt_2[1] * (a_pt_2[0] * a_pt_2[0]))) *
      a_signed_area / 60.0;
  moments.centroid()[2] =
      (3. * (a_aligned_paraboloid.a() * a_aligned_paraboloid.a()) *
           (a_pt_0[0] * a_pt_0[0] * a_pt_0[0] * a_pt_0[0] +
            a_pt_1[0] * a_pt_1[0] * a_pt_1[0] * a_pt_1[0] +
            a_pt_1[0] * a_pt_1[0] * a_pt_1[0] * a_pt_2[0] +
            a_pt_1[0] * a_pt_1[0] * (a_pt_2[0] * a_pt_2[0]) +
            a_pt_1[0] * (a_pt_2[0] * a_pt_2[0] * a_pt_2[0]) +
            a_pt_2[0] * a_pt_2[0] * a_pt_2[0] * a_pt_2[0] +
            a_pt_0[0] * a_pt_0[0] * a_pt_0[0] * (a_pt_1[0] + a_pt_2[0]) +
            a_pt_0[0] * a_pt_0[0] *
                (a_pt_1[0] * a_pt_1[0] + a_pt_1[0] * a_pt_2[0] +
                 a_pt_2[0] * a_pt_2[0]) +
            a_pt_0[0] * (a_pt_1[0] * a_pt_1[0] * a_pt_1[0] +
                         a_pt_1[0] * a_pt_1[0] * a_pt_2[0] +
                         a_pt_1[0] * (a_pt_2[0] * a_pt_2[0]) +
                         a_pt_2[0] * a_pt_2[0] * a_pt_2[0])) +
       3. * (a_aligned_paraboloid.b() * a_aligned_paraboloid.b()) *
           (a_pt_0[1] * a_pt_0[1] * a_pt_0[1] * a_pt_0[1] +
            a_pt_1[1] * a_pt_1[1] * a_pt_1[1] * a_pt_1[1] +
            a_pt_1[1] * a_pt_1[1] * a_pt_1[1] * a_pt_2[1] +
            a_pt_1[1] * a_pt_1[1] * (a_pt_2[1] * a_pt_2[1]) +
            a_pt_1[1] * (a_pt_2[1] * a_pt_2[1] * a_pt_2[1]) +
            a_pt_2[1] * a_pt_2[1] * a_pt_2[1] * a_pt_2[1] +
            a_pt_0[1] * a_pt_0[1] * a_pt_0[1] * (a_pt_1[1] + a_pt_2[1]) +
            a_pt_0[1] * a_pt_0[1] *
                (a_pt_1[1] * a_pt_1[1] + a_pt_1[1] * a_pt_2[1] +
                 a_pt_2[1] * a_pt_2[1]) +
            a_pt_0[1] * (a_pt_1[1] * a_pt_1[1] * a_pt_1[1] +
                         a_pt_1[1] * a_pt_1[1] * a_pt_2[1] +
                         a_pt_1[1] * (a_pt_2[1] * a_pt_2[1]) +
                         a_pt_2[1] * a_pt_2[1] * a_pt_2[1])) +
       a_aligned_paraboloid.a() * a_aligned_paraboloid.b() *
           (a_pt_1[0] * a_pt_2[0] *
                (a_pt_0[1] * a_pt_0[1] + 3. * (a_pt_1[1] * a_pt_1[1]) +
                 4. * a_pt_1[1] * a_pt_2[1] + 3. * (a_pt_2[1] * a_pt_2[1]) +
                 2. * a_pt_0[1] * (a_pt_1[1] + a_pt_2[1])) +
            a_pt_0[0] * a_pt_0[0] *
                (6. * (a_pt_0[1] * a_pt_0[1]) + a_pt_1[1] * a_pt_1[1] +
                 a_pt_1[1] * a_pt_2[1] + a_pt_2[1] * a_pt_2[1] +
                 3. * a_pt_0[1] * (a_pt_1[1] + a_pt_2[1])) +
            a_pt_1[0] * a_pt_1[0] *
                (a_pt_0[1] * a_pt_0[1] + 6. * (a_pt_1[1] * a_pt_1[1]) +
                 3. * a_pt_1[1] * a_pt_2[1] + a_pt_2[1] * a_pt_2[1] +
                 a_pt_0[1] * (3. * a_pt_1[1] + a_pt_2[1])) +
            a_pt_2[0] * a_pt_2[0] *
                (a_pt_0[1] * a_pt_0[1] + a_pt_1[1] * a_pt_1[1] +
                 3. * a_pt_1[1] * a_pt_2[1] + 6. * (a_pt_2[1] * a_pt_2[1]) +
                 a_pt_0[1] * (a_pt_1[1] + 3. * a_pt_2[1])) +
            a_pt_0[0] *
                (a_pt_1[0] *
                     (3. * (a_pt_0[1] * a_pt_0[1]) +
                      4. * a_pt_0[1] * a_pt_1[1] +
                      3. * (a_pt_1[1] * a_pt_1[1]) +
                      2. * a_pt_0[1] * a_pt_2[1] + 2. * a_pt_1[1] * a_pt_2[1] +
                      a_pt_2[1] * a_pt_2[1]) +
                 a_pt_2[0] *
                     (3. * (a_pt_0[1] * a_pt_0[1]) +
                      2. * a_pt_0[1] * a_pt_1[1] + a_pt_1[1] * a_pt_1[1] +
                      4. * a_pt_0[1] * a_pt_2[1] + 2. * a_pt_1[1] * a_pt_2[1] +
                      3. * (a_pt_2[1] * a_pt_2[1]))))) *
      a_signed_area / 90.0;
  return moments;
}

template <>
inline Volume computeMomentContributionClippedTriangle<Volume, 1>(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print) {
  return -(a_aligned_paraboloid.a() *
               (a_pt_0[0] * (a_pt_0[0] + a_pt_1[0] + a_pt_2[0]) +
                a_pt_1[0] * (a_pt_1[0] + a_pt_2[0]) + a_pt_2[0] * a_pt_2[0]) +
           a_aligned_paraboloid.b() *
               (a_pt_0[1] * (a_pt_0[1] + a_pt_1[1] + a_pt_2[1]) +
                a_pt_1[1] * (a_pt_1[1] + a_pt_2[1]) + a_pt_2[1] * a_pt_2[1])) *
         a_signed_area / 12.0;
}

template <>
inline Volume computeMomentContributionClippedTriangle<Volume, 2>(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print) {
  return 0.0;
}

template <>
inline Volume computeMomentContributionUnclippedTriangle<Volume, 0>(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print) {
  if (a_print) {
    amr_triangles_unclipped.insert(
        amr_triangles_unclipped.end(),
        {a_pt_0[0], a_pt_0[1], a_pt_0[2], a_pt_1[0], a_pt_1[1], a_pt_1[2],
         a_pt_2[0], a_pt_2[1], a_pt_2[2]});
  }
  return (a_pt_0[2] + a_pt_1[2] + a_pt_2[2]) * a_signed_area / 3.0;
}

template <>
inline VolumeMoments
computeMomentContributionUnclippedTriangle<VolumeMoments, 0>(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print) {
  if (a_print) {
    amr_triangles_unclipped.insert(
        amr_triangles_unclipped.end(),
        {a_pt_0[0], a_pt_0[1], a_pt_0[2], a_pt_1[0], a_pt_1[1], a_pt_1[2],
         a_pt_2[0], a_pt_2[1], a_pt_2[2]});
  }
  auto moments = VolumeMoments::fromScalarConstant(0.0);
  moments.volume() = (a_pt_0[2] + a_pt_1[2] + a_pt_2[2]) * a_signed_area / 3.0;
  moments.centroid()[0] = (2. * a_pt_0[0] * a_pt_0[2] + a_pt_1[0] * a_pt_0[2] +
                           a_pt_2[0] * a_pt_0[2] + a_pt_0[0] * a_pt_1[2] +
                           2. * a_pt_1[0] * a_pt_1[2] + a_pt_2[0] * a_pt_1[2] +
                           a_pt_0[0] * a_pt_2[2] + a_pt_1[0] * a_pt_2[2] +
                           2. * a_pt_2[0] * a_pt_2[2]) *
                          a_signed_area / 12.0;
  moments.centroid()[1] = (2. * a_pt_0[1] * a_pt_0[2] + a_pt_1[1] * a_pt_0[2] +
                           a_pt_2[1] * a_pt_0[2] + a_pt_0[1] * a_pt_1[2] +
                           2. * a_pt_1[1] * a_pt_1[2] + a_pt_2[1] * a_pt_1[2] +
                           a_pt_0[1] * a_pt_2[2] + a_pt_1[1] * a_pt_2[2] +
                           2. * a_pt_2[1] * a_pt_2[2]) *
                          a_signed_area / 12.0;
  moments.centroid()[2] =
      (a_pt_0[2] * a_pt_0[2] + a_pt_0[2] * a_pt_1[2] + a_pt_1[2] * a_pt_1[2] +
       a_pt_0[2] * a_pt_2[2] + a_pt_1[2] * a_pt_2[2] + a_pt_2[2] * a_pt_2[2]) *
      a_signed_area / 12.0;
  return moments;
}

template <>
inline Volume computeMomentContributionUnclippedTriangle<Volume, 1>(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print) {
  return (4.0 * (a_pt_0[2] + a_pt_1[2] + a_pt_2[2]) +
          (a_aligned_paraboloid.a() *
               (a_pt_0[0] * (a_pt_0[0] + a_pt_1[0] + a_pt_2[0]) +
                a_pt_1[0] * (a_pt_1[0] + a_pt_2[0]) + a_pt_2[0] * a_pt_2[0]) +
           a_aligned_paraboloid.b() *
               (a_pt_0[1] * (a_pt_0[1] + a_pt_1[1] + a_pt_2[1]) +
                a_pt_1[1] * (a_pt_1[1] + a_pt_2[1]) + a_pt_2[1] * a_pt_2[1]))) *
         a_signed_area / 12.0;
}

template <>
inline Volume computeMomentContributionUnclippedTriangle<Volume, 2>(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const double a_signed_area,
    const bool a_print) {
  return (2.0 * (a_pt_0[2] + a_pt_1[2] + a_pt_2[2]) +
          (a_aligned_paraboloid.a() *
               (a_pt_0[0] * (a_pt_0[0] + a_pt_1[0] + a_pt_2[0]) +
                a_pt_1[0] * (a_pt_1[0] + a_pt_2[0]) + a_pt_2[0] * a_pt_2[0]) +
           a_aligned_paraboloid.b() *
               (a_pt_0[1] * (a_pt_0[1] + a_pt_1[1] + a_pt_2[1]) +
                a_pt_1[1] * (a_pt_1[1] + a_pt_2[1]) + a_pt_2[1] * a_pt_2[1]))) *
         a_signed_area / 6.0;
}

template <class ReturnType>
void computeMomentContributionMixedTriangle(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const Normal& a_normal,
    const double a_signed_area,
    std::array<ReturnType, N_AMR_STRATEGIES>& a_moments_to_add,
    const bool a_print) {
  auto below =
      StackVector<bool, 3>({vertexBelow(a_pt_0, a_aligned_paraboloid),
                            vertexBelow(a_pt_1, a_aligned_paraboloid),
                            vertexBelow(a_pt_2, a_aligned_paraboloid)});

  if (below[0] && below[1] && below[2]) {
    a_moments_to_add[0] =
        computeMomentContributionUnclippedTriangle<ReturnType, 0>(
            a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
            a_print);
    if constexpr (N_AMR_STRATEGIES > 1) {
      a_moments_to_add[1] =
          computeMomentContributionUnclippedTriangle<ReturnType, 1>(
              a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
              a_print);
    }
    if constexpr (N_AMR_STRATEGIES > 2) {
      a_moments_to_add[2] =
          computeMomentContributionUnclippedTriangle<ReturnType, 2>(
              a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
              a_print);
    }
  } else if (!below[0] && !below[1] && !below[2]) {
    a_moments_to_add[0] =
        computeMomentContributionClippedTriangle<ReturnType, 0>(
            a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
            a_print);
    if constexpr (N_AMR_STRATEGIES > 1) {
      a_moments_to_add[1] =
          computeMomentContributionClippedTriangle<ReturnType, 1>(
              a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
              a_print);
    }
    if constexpr (N_AMR_STRATEGIES > 2) {
      a_moments_to_add[2] =
          computeMomentContributionClippedTriangle<ReturnType, 2>(
              a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
              a_print);
    }
  } else {
    const double d = a_normal * a_pt_0;
    auto tris =
        std::array<Pt, 6>({a_pt_0, a_pt_1, a_pt_2, a_pt_0, a_pt_1, a_pt_2});
    for (UnsignedIndex_t v = 0; v < 3; ++v) {
      tris[v][2] = -a_aligned_paraboloid.a() * tris[v][0] * tris[v][0] -
                   a_aligned_paraboloid.b() * tris[v][1] * tris[v][1];
    }

    static constexpr std::array<std::array<UnsignedIndex_t, 3>, 3> id_lists{
        {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}}};

    for (UnsignedIndex_t v = 0; v < 3; ++v) {
      const auto& ids = id_lists[v];
      if ((below[ids[0]] && !below[ids[1]] && !below[ids[2]]) ||
          (!below[ids[0]] && below[ids[1]] && below[ids[2]])) {
        const double nx0 = a_normal * tris[ids[0]];
        const double nx1 = a_normal * tris[ids[1]];
        const double nx2 = a_normal * tris[ids[2]];
        if (std::fabs(nx0 - nx1) < DBL_EPSILON ||
            std::fabs(nx0 - nx2) < DBL_EPSILON) {
          Pt center = (tris[ids[0]] + tris[ids[1]] + tris[ids[2]]);
          center /= 3.0;
          if (vertexBelow(center, a_aligned_paraboloid)) {
            a_moments_to_add[0] =
                computeMomentContributionUnclippedTriangle<ReturnType, 0>(
                    a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
                    a_print);
            if constexpr (N_AMR_STRATEGIES > 1) {
              a_moments_to_add[1] =
                  computeMomentContributionUnclippedTriangle<ReturnType, 1>(
                      a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2,
                      a_signed_area, a_print);
            }
            if constexpr (N_AMR_STRATEGIES > 2) {
              a_moments_to_add[2] =
                  computeMomentContributionUnclippedTriangle<ReturnType, 2>(
                      a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2,
                      a_signed_area, a_print);
            }
          } else {
            a_moments_to_add[0] =
                computeMomentContributionClippedTriangle<ReturnType, 0>(
                    a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
                    a_print);
            if constexpr (N_AMR_STRATEGIES > 1) {
              a_moments_to_add[1] =
                  computeMomentContributionClippedTriangle<ReturnType, 1>(
                      a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2,
                      a_signed_area, a_print);
            }
            if constexpr (N_AMR_STRATEGIES > 2) {
              a_moments_to_add[2] =
                  computeMomentContributionClippedTriangle<ReturnType, 2>(
                      a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2,
                      a_signed_area, a_print);
            }
          }
          break;
        } else {
          double t1 = (nx0 - d) / (nx0 - nx1);
          double t2 = (nx0 - d) / (nx0 - nx2);
          t1 = std::min(1.0, std::max(0.0, t1));
          t2 = std::min(1.0, std::max(0.0, t2));
          const Pt new_pt_1 = tris[ids[0]] * (1.0 - t1) + t1 * tris[ids[1]];
          const Pt new_pt_2 = tris[ids[0]] * (1.0 - t2) + t2 * tris[ids[2]];
          const double signed_area_1 =
              triangleSignedArea(tris[3 + ids[0]], new_pt_1, new_pt_2);
          const double signed_area_2 =
              triangleSignedArea(new_pt_1, tris[ids[1]], tris[ids[2]]);
          const double signed_area_3 =
              triangleSignedArea(tris[ids[2]], new_pt_2, new_pt_1);
          if (below[ids[0]]) {
            a_moments_to_add[0] =
                computeMomentContributionUnclippedTriangle<ReturnType, 0>(
                    a_aligned_paraboloid, tris[3 + ids[0]], new_pt_1, new_pt_2,
                    signed_area_1, a_print) +
                computeMomentContributionClippedTriangle<ReturnType, 0>(
                    a_aligned_paraboloid, new_pt_1, tris[3 + ids[1]],
                    tris[3 + ids[2]], signed_area_2, a_print) +
                computeMomentContributionClippedTriangle<ReturnType, 0>(
                    a_aligned_paraboloid, tris[3 + ids[2]], new_pt_2, new_pt_1,
                    signed_area_3, a_print);
            if constexpr (N_AMR_STRATEGIES > 1) {
              a_moments_to_add[1] =
                  computeMomentContributionUnclippedTriangle<ReturnType, 1>(
                      a_aligned_paraboloid, tris[3 + ids[0]], new_pt_1,
                      new_pt_2, signed_area_1, a_print) +
                  computeMomentContributionClippedTriangle<ReturnType, 1>(
                      a_aligned_paraboloid, new_pt_1, tris[3 + ids[1]],
                      tris[3 + ids[2]], signed_area_2, a_print) +
                  computeMomentContributionClippedTriangle<ReturnType, 1>(
                      a_aligned_paraboloid, tris[3 + ids[2]], new_pt_2,
                      new_pt_1, signed_area_3, a_print);
            }
            if constexpr (N_AMR_STRATEGIES > 2) {
              a_moments_to_add[2] =
                  computeMomentContributionUnclippedTriangle<ReturnType, 2>(
                      a_aligned_paraboloid, tris[3 + ids[0]], new_pt_1,
                      new_pt_2, signed_area_1, a_print) +
                  computeMomentContributionClippedTriangle<ReturnType, 2>(
                      a_aligned_paraboloid, new_pt_1, tris[3 + ids[1]],
                      tris[3 + ids[2]], signed_area_2, a_print) +
                  computeMomentContributionClippedTriangle<ReturnType, 2>(
                      a_aligned_paraboloid, tris[3 + ids[2]], new_pt_2,
                      new_pt_1, signed_area_3, a_print);
            }
          } else {
            a_moments_to_add[0] =
                computeMomentContributionClippedTriangle<ReturnType, 0>(
                    a_aligned_paraboloid, tris[3 + ids[0]], new_pt_1, new_pt_2,
                    signed_area_1, a_print) +
                computeMomentContributionUnclippedTriangle<ReturnType, 0>(
                    a_aligned_paraboloid, new_pt_1, tris[3 + ids[1]],
                    tris[3 + ids[2]], signed_area_2, a_print) +
                computeMomentContributionUnclippedTriangle<ReturnType, 0>(
                    a_aligned_paraboloid, tris[3 + ids[2]], new_pt_2, new_pt_1,
                    signed_area_3, a_print);
            if constexpr (N_AMR_STRATEGIES > 1) {
              a_moments_to_add[1] =
                  computeMomentContributionClippedTriangle<ReturnType, 1>(
                      a_aligned_paraboloid, tris[3 + ids[0]], new_pt_1,
                      new_pt_2, signed_area_1, a_print) +
                  computeMomentContributionUnclippedTriangle<ReturnType, 1>(
                      a_aligned_paraboloid, new_pt_1, tris[3 + ids[1]],
                      tris[3 + ids[2]], signed_area_2, a_print) +
                  computeMomentContributionUnclippedTriangle<ReturnType, 1>(
                      a_aligned_paraboloid, tris[3 + ids[2]], new_pt_2,
                      new_pt_1, signed_area_3, a_print);
            }
            if constexpr (N_AMR_STRATEGIES > 2) {
              a_moments_to_add[2] =
                  computeMomentContributionClippedTriangle<ReturnType, 2>(
                      a_aligned_paraboloid, tris[3 + ids[0]], new_pt_1,
                      new_pt_2, signed_area_1, a_print) +
                  computeMomentContributionUnclippedTriangle<ReturnType, 2>(
                      a_aligned_paraboloid, new_pt_1, tris[3 + ids[1]],
                      tris[3 + ids[2]], signed_area_2, a_print) +
                  computeMomentContributionUnclippedTriangle<ReturnType, 2>(
                      a_aligned_paraboloid, tris[3 + ids[2]], new_pt_2,
                      new_pt_1, signed_area_3, a_print);
            }
          }
          break;
        }
      }
      assert(v < 2);
    }
  }
}

std::pair<bool, bool> computeZBounds(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2) {
  // Compute bounding box of triangle
  std::array<double, 6> tri_bounds;
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    tri_bounds[2 * d] = std::min({a_pt_0[d], a_pt_1[d], a_pt_2[d]});
    tri_bounds[2 * d + 1] = std::max({a_pt_0[d], a_pt_1[d], a_pt_2[d]});
  }
  //   if (tri_bounds[0] * tri_bounds[1] < 0.0 ||
  //       tri_bounds[2] * tri_bounds[3] < 0.0) {
  // Compute max gradient magnitude of paraboloid over triangle bounding box
  const double max_grad_x =
      std::max(std::fabs(2.0 * a_aligned_paraboloid.a() * tri_bounds[0]),
               std::fabs(2.0 * a_aligned_paraboloid.a() * tri_bounds[1]));
  const double max_grad_y =
      std::max(std::fabs(2.0 * a_aligned_paraboloid.b() * tri_bounds[2]),
               std::fabs(2.0 * a_aligned_paraboloid.b() * tri_bounds[3]));
  // Compute paraboloid height in middle of box
  const double z_paraboloid =
      -0.25 * (a_aligned_paraboloid.a() * (tri_bounds[0] + tri_bounds[1]) *
                   (tri_bounds[0] + tri_bounds[1]) +
               a_aligned_paraboloid.b() * (tri_bounds[2] + tri_bounds[3]) *
                   (tri_bounds[2] + tri_bounds[3]));
  // Compute bouding box of paraboloid (over-estimate)
  const double para_max =
      z_paraboloid + 0.5 * (max_grad_x * (tri_bounds[1] - tri_bounds[0]) +
                            max_grad_y * (tri_bounds[3] - tri_bounds[2]));
  const double para_min = 2.0 * z_paraboloid - para_max;
  // Return min/max of z height of triangle and paraboloid
  return std::pair<bool, bool>({para_min > tri_bounds[5] + AMR_DBL_EPSILON,
                                tri_bounds[4] > para_max + AMR_DBL_EPSILON});
  //   } else {
  //     const double z00 =
  //         -(a_aligned_paraboloid.a() * tri_bounds[0] * tri_bounds[0] +
  //           a_aligned_paraboloid.b() * tri_bounds[2] * tri_bounds[2]);
  //     const double z01 =
  //         -(a_aligned_paraboloid.a() * tri_bounds[0] * tri_bounds[0] +
  //           a_aligned_paraboloid.b() * tri_bounds[3] * tri_bounds[3]);
  //     const double z10 =
  //         -(a_aligned_paraboloid.a() * tri_bounds[1] * tri_bounds[1] +
  //           a_aligned_paraboloid.b() * tri_bounds[2] * tri_bounds[2]);
  //     const double z11 =
  //         -(a_aligned_paraboloid.a() * tri_bounds[1] * tri_bounds[1] +
  //           a_aligned_paraboloid.b() * tri_bounds[3] * tri_bounds[3]);
  //     const double para_max = std::max({z00, z01, z10, z11});
  //     const double para_min = std::min({z00, z01, z10, z11});
  //     // Return min/max of z height of triangle and paraboloid
  //     return std::pair<bool, bool>({para_min > tri_bounds[5] +
  //     AMR_DBL_EPSILON,
  //                                   tri_bounds[4] > para_max +
  //                                   AMR_DBL_EPSILON});
  //   }
}

template <class ReturnType>
void computeMomentContributionAMR(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1, const Pt& a_pt_2, const Normal& a_normal,
    const double a_signed_area, const UnsignedIndex_t a_amr_level,
    const UnsignedIndex_t a_max_amr_level,
    std::array<std::pair<ReturnType, ReturnType>, N_AMR_STRATEGIES>&
        a_full_moments_ref,
    std::array<std::pair<ReturnType, ReturnType>, N_AMR_STRATEGIES>&
        a_full_moments,
    const bool a_print) {
  std::array<ReturnType, N_AMR_STRATEGIES> moments_to_add;

  // Compute z-bounds of triangle and paraboloid
  auto z_limits = computeZBounds(a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2);
  if (z_limits.first) {
    // Max of triangle is smaller than min of paraboloid
    moments_to_add[0] =
        computeMomentContributionUnclippedTriangle<ReturnType, 0>(
            a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
            a_print);

    if constexpr (N_AMR_STRATEGIES > 1) {
      moments_to_add[1] =
          computeMomentContributionUnclippedTriangle<ReturnType, 1>(
              a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
              a_print);
    }

    if constexpr (N_AMR_STRATEGIES > 2) {
      moments_to_add[2] =
          computeMomentContributionUnclippedTriangle<ReturnType, 2>(
              a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
              a_print);
    }
    kahanSummationMoments<ReturnType>(a_full_moments, a_full_moments_ref,
                                      moments_to_add);
    return;
  } else if (z_limits.second) {
    // Max of triangle is smaller than min of paraboloid
    moments_to_add[0] = computeMomentContributionClippedTriangle<ReturnType, 0>(
        a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area, a_print);

    if constexpr (N_AMR_STRATEGIES > 1) {
      moments_to_add[1] =
          computeMomentContributionClippedTriangle<ReturnType, 1>(
              a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
              a_print);
    }

    if constexpr (N_AMR_STRATEGIES > 2) {
      moments_to_add[2] =
          computeMomentContributionClippedTriangle<ReturnType, 2>(
              a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_signed_area,
              a_print);
    }
    kahanSummationMoments<ReturnType>(a_full_moments, a_full_moments_ref,
                                      moments_to_add);
    return;
  } else if (a_amr_level == a_max_amr_level) {
    computeMomentContributionMixedTriangle<ReturnType>(
        a_aligned_paraboloid, a_pt_0, a_pt_1, a_pt_2, a_normal, a_signed_area,
        moments_to_add, a_print);

    kahanSummationMoments<ReturnType>(a_full_moments, a_full_moments_ref,
                                      moments_to_add);
    return;
  }
  // Refine triangle and call function recursively
  const Pt c0 = 0.5 * (a_pt_0 + a_pt_1);
  const Pt c1 = 0.5 * (a_pt_1 + a_pt_2);
  const Pt c2 = 0.5 * (a_pt_2 + a_pt_0);
  computeMomentContributionAMR<ReturnType>(
      a_aligned_paraboloid, a_pt_0, c0, c2, a_normal, 0.25 * a_signed_area,
      a_amr_level + 1, a_max_amr_level, a_full_moments_ref, a_full_moments,
      a_print);
  computeMomentContributionAMR<ReturnType>(
      a_aligned_paraboloid, a_pt_1, c1, c0, a_normal, 0.25 * a_signed_area,
      a_amr_level + 1, a_max_amr_level, a_full_moments_ref, a_full_moments,
      a_print);
  computeMomentContributionAMR<ReturnType>(
      a_aligned_paraboloid, a_pt_2, c2, c1, a_normal, 0.25 * a_signed_area,
      a_amr_level + 1, a_max_amr_level, a_full_moments_ref, a_full_moments,
      a_print);
  computeMomentContributionAMR<ReturnType>(
      a_aligned_paraboloid, c0, c1, c2, a_normal, 0.25 * a_signed_area,
      a_amr_level + 1, a_max_amr_level, a_full_moments_ref, a_full_moments,
      a_print);
}

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithParaboloidAMR(
    SegmentedHalfEdgePolyhedronType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const AlignedParaboloid& a_paraboloid,
    const UnsignedIndex_t a_max_amr_level, const std::string& a_filename) {
  const bool print = !a_filename.empty();
  // We have N strategies for computing the moments (and will later on choose
  // the best of those)
  std::array<std::pair<ReturnType, ReturnType>, N_AMR_STRATEGIES> full_moments;
  std::array<std::pair<ReturnType, ReturnType>, N_AMR_STRATEGIES>
      full_moments_ref;
  for (UnsignedIndex_t m = 0; m < N_AMR_STRATEGIES; ++m) {
    full_moments[m] = std::pair<ReturnType, ReturnType>(
        {ReturnType::fromScalarConstant(0.0),
         ReturnType::fromScalarConstant(0.0)});
    full_moments_ref[m] = std::pair<ReturnType, ReturnType>(
        {ReturnType::fromScalarConstant(0.0),
         ReturnType::fromScalarConstant(0.0)});
  }
  // We loop over each face, triangulate, and split until there are no
  // intersections OR we have reached the max refinement level
  const auto number_of_faces = a_polytope->getNumberOfFaces();
  for (UnsignedIndex_t f = 0; f < number_of_faces; ++f) {
    const auto& face = *(*a_polytope)[f];
    const auto& face_normal = face.getPlane().normal();
    const auto starting_half_edge = face.getStartingHalfEdge();
    const auto& ref_pt = starting_half_edge->getVertex()->getLocation().getPt();
    auto current_half_edge =
        starting_half_edge->getNextHalfEdge()->getNextHalfEdge();
    auto prev_pt =
        current_half_edge->getPreviousVertex()->getLocation().getPt();
    do {
      const auto& curr_pt =
          current_half_edge->getVertex()->getLocation().getPt();
      const double signed_area = triangleSignedArea(ref_pt, prev_pt, curr_pt);
      computeMomentContributionAMR<ReturnType>(
          a_paraboloid, ref_pt, prev_pt, curr_pt, face_normal, signed_area, 0,
          a_max_amr_level, full_moments_ref, full_moments, print);
      prev_pt = curr_pt;
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != starting_half_edge);
  }

  // Print triangles
  if (print) {
    auto write_triangles = [a_filename](
                               const std::string& a_header,
                               const std::string& a_file_prefix,
                               const std::vector<double>& a_triangle_list) {
      // binary file
      char head[80];
      std::strncpy(head, a_header.c_str(), a_header.size() - 1);
      char attribute[2] = "0";
      char dummy[4] = "0";
      const auto nTriLong = amr_triangles_clipped.size() / 3;

      std::ofstream myfile;

      myfile.open(a_file_prefix + a_filename + ".stl",
                  std::ios::out | std::ios::binary);
      myfile.write(head, sizeof(head));
      myfile.write(reinterpret_cast<const char*>(&nTriLong), 4);

      // write down every triangle
      for (std::size_t i = 0; i < a_triangle_list.size(); i += 9) {
        // normal vector coordinates
        myfile.write(dummy, 4);
        myfile.write(dummy, 4);
        myfile.write(dummy, 4);
        for (std::size_t t = 0; t < 9; ++t) {
          const float val = static_cast<float>(a_triangle_list[i + t]);
          myfile.write(reinterpret_cast<const char*>(&val), 4);
        }
        myfile.write(attribute, 2);
      }
    };
    write_triangles("clipped_triangles", "clipped", amr_triangles_clipped);
    write_triangles("unclipped_triangles", "unclipped",
                    amr_triangles_unclipped);
  }

  // Find best estimate
  auto best_full_moments = full_moments[0];
  for (UnsignedIndex_t m = 1; m < full_moments.size() - 1; ++m) {
    if constexpr (std::is_same<ReturnType, Volume>::value) {
      if (full_moments[m].second < best_full_moments.second) {
        best_full_moments = full_moments[m];
      }
    } else if constexpr (std::is_same<ReturnType, VolumeMoments>::value) {
      if (full_moments[m].second.volume() < best_full_moments.second.volume()) {
        best_full_moments = full_moments[m];
      }
    }
  }

  amr_triangles_unclipped.clear();
  amr_triangles_clipped.clear();

  return best_full_moments.first;
}
}  // namespace IRL
#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_AMR_TPP_
