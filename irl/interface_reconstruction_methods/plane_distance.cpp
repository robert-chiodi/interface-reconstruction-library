// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/interface_reconstruction_methods/plane_distance.h"

#include "irl/helpers/helper.h"
#include "irl/parameters/constants.h"

namespace IRL {

namespace plane_distance_details {

static double getAlphaRectangularCuboid(const std::array<double, 3>& a_mm,
                                        const double a_VOFo) {
  // Copy and take absolute value of mm values
  auto m = std::array<double, 3>{std::fabs(a_mm[0]), std::fabs(a_mm[1]),
                                 std::fabs(a_mm[2])};
  sort3Ascending(m.data());
  const double m12 = {m[0] + m[1]};
  const double V1 = {m[0] * m[0] / safelyEpsilon(6.0 * m[1] * m[2])};
  const double V2 = {V1 + 0.5 * (m[1] - m[0]) / m[2]};
  const double V3 = m12 <= m[2] ? 0.5 * m12 / m[2]
                                : (m[2] * m[2] * (3.0 * m12 - m[2]) +
                                   m[0] * m[0] * (m[0] - 3.0 * m[2]) +
                                   m[1] * m[1] * (m[1] - 3.0 * m[2])) /
                                      safelyEpsilon(6.0 * m[0] * m[1] * m[2]);
  if (a_VOFo >= 0.0 && a_VOFo < V1) {
    return std::pow((6.0 * m[0] * m[1] * m[2] * a_VOFo), 1.0 / 3.0);
  } else if (a_VOFo >= V1 && a_VOFo < V2) {
    return 0.5 * (m[0] + sqrt(m[0] * m[0] + 8.0 * m[1] * m[2] * (a_VOFo - V1)));
  } else if (a_VOFo >= V2 && a_VOFo < V3) {
    double a0 = -(m[0] * m[0] * m[0] + m[1] * m[1] * m[1] -
                  6.0 * m[0] * m[1] * m[2] * a_VOFo);
    double a1 = 3.0 * (m[0] * m[0] + m[1] * m[1]);
    double a2 = -3.0 * m12;
    double p0 = -(a1 / 3.0 - a2 * a2 / 9.0);
    double q0 = (a1 * a2 - 3.0 * a0) / 6.0 - a2 * a2 * a2 / 27.0;
    double theta =
        std::acos(clipBetween(-1.0, q0 / std::sqrt(p0 * p0 * p0), 1.0)) / 3.0;
    return std::sqrt(p0) *
               (std::sqrt(3.0) * std::sin(theta) - std::cos(theta)) -
           a2 / 3.0;
  } else {
    if (m[2] >= m12) {
      return m[2] * a_VOFo + 0.5 * m12;
    } else {
      double a0 =
          -0.5 * (m[0] * m[0] * m[0] + m[1] * m[1] * m[1] + m[2] * m[2] * m[2] -
                  6.0 * m[0] * m[1] * m[2] * a_VOFo);
      double a1 = 1.5 * (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
      double a2 = -1.5;
      double p0 = -(a1 / 3.0 - a2 * a2 / 9.0);
      double q0 = (a1 * a2 - 3.0 * a0) / 6.0 - a2 * a2 * a2 / 27.0;
      double theta =
          std::acos(clipBetween(-1.0, q0 / std::sqrt(p0 * p0 * p0), 1.0)) / 3.0;
      return std::sqrt(p0) *
                 (std::sqrt(3.0) * std::sin(theta) - std::cos(theta)) -
             a2 / 3.0;
    }
  }
}

}  // namespace plane_distance_details

double findDistanceOnePlane(const RectangularCuboid& a_rectangular_cuboid,
                            const double a_volume_fraction,
                            const Normal& a_normal) {
  auto mm = std::array<double, 3>{
      {a_normal[0] * a_rectangular_cuboid.calculateSideLength(0),
       a_normal[1] * a_rectangular_cuboid.calculateSideLength(1),
       a_normal[2] * a_rectangular_cuboid.calculateSideLength(2)}};
  double norm = std::fabs(mm[0]) + std::fabs(mm[1]) + std::fabs(mm[2]);
  double factor = {0.0};
  for (auto& elem : mm) {
    elem /= norm;
    if (elem < 0.0) {
      factor += elem;
    }
  }
  double VOFo =
      a_volume_fraction > 0.5 ? 1.0 - a_volume_fraction : a_volume_fraction;
  double alpha = plane_distance_details::getAlphaRectangularCuboid(mm, VOFo);
  alpha = a_volume_fraction > 0.5 ? 1.0 - alpha : alpha;
  alpha += (factor - 0.5 * (mm[0] + mm[1] + mm[2]));
  return alpha * norm + a_normal * a_rectangular_cuboid.calculateCentroid();
}

namespace plane_distance_details {
static UnsignedIndex_t getRegime(const double a_volume_fraction,
                                 const std::array<double, 4>& a_vd,
                                 double* a_section_fraction) {
  if (a_vd[1] < DBL_MIN && a_vd[2] < DBL_MIN) {
    *a_section_fraction = 1.0;
    return 2;
  }

  if (a_vd[2] - a_vd[1] < DBL_MIN && a_vd[3] - a_vd[2] < DBL_MIN) {
    *a_section_fraction = 1.0;
    return 0;
  }

  // Check if in regime 0
  const double volume_fraction_b = (a_vd[1] / a_vd[2]) * (a_vd[1] / a_vd[3]);
  if (a_volume_fraction <= volume_fraction_b) {
    *a_section_fraction = volume_fraction_b;
    return 0;
  }

  // Check if in regime 2
  const double volume_fraction_c =
      (a_vd[3] - a_vd[2]) * (a_vd[3] - a_vd[2]) / (a_vd[3] - a_vd[1]) / a_vd[3];
  if (a_volume_fraction >= (1.0 - volume_fraction_c)) {
    *a_section_fraction = volume_fraction_c;
    return 2;
  } else {
    return 1;
  }

  return static_cast<UnsignedIndex_t>(-1);  // Should never reach this.
}

template <UnsignedIndex_t kRegimeNumber>
static double calculateRegimeDistance(const Tet& a_tet,
                                      const double a_volume_fraction,
                                      const Normal& a_normal,
                                      const double a_section_volume_fraction,
                                      const std::array<double, 4>& a_vd);

template <>
double calculateRegimeDistance<0>(const Tet& a_tet,
                                  const double a_volume_fraction,
                                  const Normal& a_normal,
                                  const double a_section_volume_fraction,
                                  const std::array<double, 4>& a_vd) {
  return a_vd[0] +
         std::pow(a_volume_fraction / a_section_volume_fraction, 1.0 / 3.0) *
             (a_vd[1] - a_vd[0]);
}

template <>
double calculateRegimeDistance<1>(const Tet& a_tet,
                                  const double a_volume_fraction,
                                  const Normal& a_normal,
                                  const double a_section_volume_fraction,
                                  const std::array<double, 4>& a_vd) {
  const double pc_m_pb = a_vd[2] - a_vd[1];
  const double pd_m_pb = a_vd[3] - a_vd[1];
  const double pd_t_pc = a_vd[3] * a_vd[2];

  const double a =
      -pc_m_pb * pc_m_pb / a_vd[3] * (1.0 / a_vd[2] + 1.0 / pd_m_pb);

  const double b = 3.0 * pc_m_pb * pc_m_pb / pd_t_pc;

  const double c = 3.0 * a_vd[1] * pc_m_pb / pd_t_pc;

  const double d = a_vd[1] * a_vd[1] / pd_t_pc - a_volume_fraction;

  const double p = c / a - b * b / (3.0 * a * a);
  const double q =
      d / a + 2.0 * b * b * b / (27.0 * a * a * a) - b * c / (3.0 * a * a);
  const double theta =
      acos(-0.5 * q / std::sqrt(-(p / 3.0 * p / 3.0 * p / 3.0))) / 3.0;

  const double alpha =
      std::sqrt(-p / 3.0) * (std::sqrt(3.0) * sin(theta) - cos(theta)) -
      b / (3.0 * a);

  return a_vd[1] + alpha * (a_vd[2] - a_vd[1]);
}

template <>
double calculateRegimeDistance<2>(const Tet& a_tet,
                                  const double a_volume_fraction,
                                  const Normal& a_normal,
                                  const double a_section_volume_fraction,
                                  const std::array<double, 4>& a_vd) {
  return a_vd[3] +
         std::pow((1.0 - a_volume_fraction) / a_section_volume_fraction,
                  1.0 / 3.0) *
             (a_vd[2] - a_vd[3]);
}

}  // namespace plane_distance_details

// Analytical distance finding on a Tet from
// Yang & James. "Analytic relations for reconstruction piecewise linear
// interfaces in triangular and tetrahedral grids. Journal of Computational
// Physics, 2006.
//
// Note, in the implementation, most things are converted from
// one-indexing to zero-indexing. Also, only mapping to B/C/D vertices
// and their distances are done, since everything is calculated relative to
// vertex A, so that will always be a_tet[0].
double findDistanceOnePlane(const Tet& a_tet, const double a_volume_fraction,
                            const Normal& a_normal) {
  std::array<double, 4> vertex_dist{{a_normal * a_tet[0], a_normal * a_tet[1],
                                     a_normal * a_tet[2], a_normal * a_tet[3]}};

  std::sort(vertex_dist.begin(), vertex_dist.end());
  std::array<double, 4> sorted_distance_to_plane;
  const double lowest_value = vertex_dist[0];
  for (auto& elem : vertex_dist) {
    elem -= lowest_value;
  }
  double section_volume_fraction;
  const auto regime = plane_distance_details::getRegime(
      a_volume_fraction, vertex_dist, &section_volume_fraction);

  switch (regime) {
    case 0:
      return plane_distance_details::calculateRegimeDistance<0>(
                 a_tet, a_volume_fraction, a_normal, section_volume_fraction,
                 vertex_dist) +
             lowest_value;
      break;

    case 1:
      return plane_distance_details::calculateRegimeDistance<1>(
                 a_tet, a_volume_fraction, a_normal, section_volume_fraction,
                 vertex_dist) +
             lowest_value;
      break;

    case 2:
      return plane_distance_details::calculateRegimeDistance<2>(
                 a_tet, a_volume_fraction, a_normal, section_volume_fraction,
                 vertex_dist) +
             lowest_value;
      break;
    default:
      std::cout
          << "Unknown regime found in Tet distance finding. Unknown regime is "
          << regime << std::endl;
      std::exit(-1);
  }
}

}  // namespace IRL
