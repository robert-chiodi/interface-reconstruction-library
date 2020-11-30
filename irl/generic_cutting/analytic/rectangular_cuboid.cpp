// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/generic_cutting/analytic/rectangular_cuboid.h"

namespace IRL {

namespace analytic_rectangular_cuboid {
static double cubeIfPositive(const double a_val) {
  return a_val > 0.0 ? a_val * a_val * a_val : 0.0;
}

static double squareIfPositive(const double a_val) {
  return a_val > 0.0 ? a_val * a_val : 0.0;
}

static double calculateVolumeByRegime(const std::array<double, 3>& a_mm,
                                      const double a_alpha) {
  assert(a_mm[0] >= 0.0);
  assert(a_mm[1] >= 0.0);
  assert(a_mm[2] >= 0.0);
  assert(a_mm[0] <= 1.0);
  assert(a_mm[1] <= 1.0);
  assert(a_mm[2] <= 1.0);
  assert(a_alpha >= 0.0);
  assert(a_alpha <= 0.5);

  const double m12 = {a_mm[0] + a_mm[1]};

  if (a_alpha >= 0.0 && a_alpha < a_mm[0]) {
    return a_alpha * a_alpha * a_alpha / (6.0 * a_mm[0] * a_mm[1] * a_mm[2]);
  } else if (a_alpha >= a_mm[0] && a_alpha < a_mm[1]) {
    const double V1 = {a_mm[0] * a_mm[0] / safelyTiny(6.0 * a_mm[1] * a_mm[2])};
    return a_alpha * (a_alpha - a_mm[0]) / (2.0 * a_mm[1] * a_mm[2]) + V1;
  } else if (a_alpha >= a_mm[1] && a_alpha < std::min(m12, a_mm[2])) {
    return (a_alpha * a_alpha * (3.0 * m12 - a_alpha) +
            a_mm[0] * a_mm[0] * (a_mm[0] - 3.0 * a_alpha) +
            a_mm[1] * a_mm[1] * (a_mm[1] - 3.0 * a_alpha)) /
           (6.0 * a_mm[0] * a_mm[1] * a_mm[2]);
  } else {
    if (a_mm[2] < m12) {
      return (a_alpha * a_alpha * (3.0 - 2.0 * a_alpha) +
              a_mm[0] * a_mm[0] * (a_mm[0] - 3.0 * a_alpha) +
              a_mm[1] * a_mm[1] * (a_mm[1] - 3.0 * a_alpha) +
              a_mm[2] * a_mm[2] * (a_mm[2] - 3.0 * a_alpha)) /
             (6.0 * a_mm[0] * a_mm[1] * a_mm[2]);
    } else {
      return (2.0 * a_alpha - m12) / (2.0 * a_mm[2]);
    }
  }
}

}  // namespace analytic_rectangular_cuboid

Volume getAnalyticVolume(const RectangularCuboid& a_rectangular_cuboid,
                         const Plane& a_plane) {
  const auto plane_normal = a_plane.normal();
  const Volume cell_volume = a_rectangular_cuboid.calculateVolume();
  if (squaredMagnitude(plane_normal) < DBL_MIN) {
    return a_plane.distance() > 0.0 ? cell_volume : Volume(0.0);
  }
  auto mm = std::array<double, 3>{
      {plane_normal[0] * a_rectangular_cuboid.calculateSideLength(0),
       plane_normal[1] * a_rectangular_cuboid.calculateSideLength(1),
       plane_normal[2] * a_rectangular_cuboid.calculateSideLength(2)}};
  const Pt lower_point = a_rectangular_cuboid.calculateCentroid() -
                         0.5 * Pt(a_rectangular_cuboid.calculateSideLength(0),
                                  a_rectangular_cuboid.calculateSideLength(1),
                                  a_rectangular_cuboid.calculateSideLength(2));
  const double norm = std::fabs(mm[0]) + std::fabs(mm[1]) + std::fabs(mm[2]);
  double alpha = (a_plane.distance() - plane_normal * lower_point) / norm;
  for (auto& element : mm) {
    element /= norm;
    if (element < 0.0) {
      element = -element;
      alpha += element;
    }
  }
  if (alpha <= 0.0) {
    return 0.0;
  } else if (alpha >= 1.0) {
    return cell_volume;
  }
  sort3Ascending(mm.data());
  const double alpha_0 = alpha > 0.5 ? 1.0 - alpha : alpha;
  const double tmp_vol =
      analytic_rectangular_cuboid::calculateVolumeByRegime(mm, alpha_0);
  return alpha > 0.5 ? (1.0 - tmp_vol) * cell_volume : tmp_vol * cell_volume;
}

}  // namespace IRL
