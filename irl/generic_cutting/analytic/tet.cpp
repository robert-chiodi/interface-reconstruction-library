// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cmath>

#include "irl/generic_cutting/analytic/tet.h"

namespace IRL {

namespace analytic_tet {
// Gets regime for the analytic tet volume calculation.
// With a_vd sorted in ascending order, regime is between the
// two vertices the intersection lays. This occurs when switching from
// negative to positive. This leads to regimes of
// Regime -1 : No intersection, Tet completely above plane
// Regime 0 : Intersectiom between X_A and X_B
// Regime 1 : Intersectiom between X_B and X_C
// Regime 2 : Intersectiom between X_C and X_D
// Regime 3 : No intersection, Tet completely below plane
static int getRegime(const std::array<double, 4>& a_vd,
                     const double a_lowest_val, double* a_section_fraction) {
  // Special cases
  if (a_vd[1] < DBL_MIN && a_vd[2] < DBL_MIN) {
    *a_section_fraction = 1.0;
    return 2;
  }

  if (a_vd[2] - a_vd[1] < DBL_MIN && a_vd[3] - a_vd[2] < DBL_MIN) {
    *a_section_fraction = 1.0;
    return 0;
  }

  // Look for first positive number to find regime
  UnsignedIndex_t index;
  for (index = 0; index < 4; ++index) {
    if (a_vd[index] + a_lowest_val >= 0.0) {
      break;
    }
  }
  const int regime = static_cast<int>(index) - 1;

  if (regime == 0) {
    *a_section_fraction = (a_vd[1] / a_vd[2]) * (a_vd[1] / a_vd[3]);

  } else if (regime == 2) {
    *a_section_fraction = (a_vd[3] - a_vd[2]) * (a_vd[3] - a_vd[2]) /
                          (a_vd[3] - a_vd[1]) / a_vd[3];
  }

  return regime;
}

template <UnsignedIndex_t kRegimeNumber>
static double calculateRegimeVolume(const std::array<double, 4>& a_vd,
                                    const double a_lowest_val,
                                    const double a_section_volume_fraction);

template <>
double calculateRegimeVolume<0>(const std::array<double, 4>& a_vd,
                                const double a_lowest_val,
                                const double a_section_volume_fraction) {
  const double mu = -(a_vd[0] + a_lowest_val) / safelyTiny(a_vd[1] - a_vd[0]);
  return std::pow(mu, 3.0) * a_section_volume_fraction;
}

template <>
double calculateRegimeVolume<1>(const std::array<double, 4>& a_vd,
                                const double a_lowest_val,
                                const double a_section_volume_fraction) {
  const double pc_m_pb = a_vd[2] - a_vd[1];
  const double pd_m_pb = a_vd[3] - a_vd[1];
  const double pd_t_pc = a_vd[3] * a_vd[2];

  const double a =
      -pc_m_pb * pc_m_pb / a_vd[3] * (1.0 / a_vd[2] + 1.0 / pd_m_pb);

  const double b = 3.0 * pc_m_pb * pc_m_pb / pd_t_pc;

  const double c = 3.0 * a_vd[1] * pc_m_pb / pd_t_pc;

  const double mu = -(a_vd[1] + a_lowest_val) / safelyTiny(a_vd[2] - a_vd[1]);
  return a_vd[1] * a_vd[1] / (a_vd[3] * a_vd[2]) + a * std::pow(mu, 3.0) +
         b * std::pow(mu, 2.0) + c * mu;
}

template <>
double calculateRegimeVolume<2>(const std::array<double, 4>& a_vd,
                                const double a_lowest_val,
                                const double a_section_volume_fraction) {
  const double mu = -(a_vd[3] + a_lowest_val) / safelyTiny(a_vd[2] - a_vd[3]);
  return std::pow(mu, 3.0) * a_section_volume_fraction;
}

}  // namespace analytic_tet

Volume getAnalyticVolume(const Tet& a_tet, const Plane& a_plane) {
  const auto& plane_normal = a_plane.normal();
if (squaredMagnitude(a_plane.normal()) < DBL_MIN) {
  return a_plane.distance() > 0.0 ? a_tet.calculateVolume() : static_cast<Volume>(0.0);
  }
  std::array<double, 4> vertex_dist{
      {plane_normal * a_tet[0] - a_plane.distance(),
       plane_normal * a_tet[1] - a_plane.distance(),
       plane_normal * a_tet[2] - a_plane.distance(),
       plane_normal * a_tet[3] - a_plane.distance()}};

  std::sort(vertex_dist.begin(), vertex_dist.end());

  const double lowest_value = vertex_dist[0];
  for (auto& elem : vertex_dist) {
    elem -= lowest_value;
  }
  double section_volume_fraction;
  const auto regime = analytic_tet::getRegime(vertex_dist, lowest_value,
                                              &section_volume_fraction);
  switch (regime) {
    case -1:
      return 0.0;
    case 0:
      return analytic_tet::calculateRegimeVolume<0>(vertex_dist, lowest_value,
                                                    section_volume_fraction) *
             a_tet.calculateVolume();
      break;

    case 1:
      return analytic_tet::calculateRegimeVolume<1>(vertex_dist, lowest_value,
                                                    section_volume_fraction) *
             a_tet.calculateVolume();
      break;

    case 2:
      return a_tet.calculateVolume() *
             (1.0 - analytic_tet::calculateRegimeVolume<2>(
                        vertex_dist, lowest_value, section_volume_fraction));
      break;
    case 3:
      return a_tet.calculateVolume();
    default:
      std::cout << "Unknown regime found in Tet distance finding. Unknown "
                   "regime is "
                << regime << std::endl;
      std::exit(-1);
  }
}

}  // namespace IRL
