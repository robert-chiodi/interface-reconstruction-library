// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/generic_cutting/analytic/rectangular_cuboid.h"
#include "irl/generic_cutting/analytic/tet.h"

#include <random>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/geometry/polyhedrons/tet.h"
#include "irl/interface_reconstruction_methods/plane_distance.h"
#include "irl/moments/volume.h"

#include "gtest/gtest.h"

namespace {

using namespace IRL;

TEST(AnalyticCutting, RectangularCuboidForVolume) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(
      global_constants::VF_LOW + DBL_EPSILON,
      global_constants::VF_HIGH - DBL_EPSILON);

  static const int ncycles = 10000;

  RectangularCuboid cuboid = unit_cell;
  cuboid =
      RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(4.0, 16.0, 2.5));
  for (auto& vertex : cuboid) {
    vertex += Pt(-10.0, 5.0, 1.5);
  }
  const auto cell_volume = cuboid.calculateVolume();

  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal correct_normal = Normal::normalized(
        random_normal(eng), random_normal(eng), random_normal(eng));
    double set_VF = random_VF(eng);
    PlanarSeparator correct_reconstruction = PlanarSeparator::fromOnePlane(
        Plane(correct_normal,
              findDistanceOnePlane(cuboid, set_VF, correct_normal)));
    auto calculated_volume =
        getAnalyticVolume(cuboid, correct_reconstruction[0]);
    EXPECT_NEAR(calculated_volume / cell_volume, set_VF, 1.0e-14)
        << correct_reconstruction << " " << set_VF << std::endl;
  }
}

TEST(AnalyticCutting, TetForVolume) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(
      global_constants::VF_LOW + DBL_EPSILON,
      global_constants::VF_HIGH - DBL_EPSILON);

  static const int ncycles = 10000;

  Tet tet({Pt(1.0, 0.0, -0.5), Pt(1.0, 1.0, 0.0), Pt(1.0, 0.0, 0.5),
           Pt(0.5, 0.0, 0.0)});
  for (auto& vertex : tet) {
    vertex += Pt(-10.0, 5.0, 1.5);
  }
  const auto tet_volume = tet.calculateVolume();

  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal correct_normal = Normal::normalized(
        random_normal(eng), random_normal(eng), random_normal(eng));
    double set_VF = random_VF(eng);
    PlanarSeparator correct_reconstruction =
        PlanarSeparator::fromOnePlane(Plane(
            correct_normal, findDistanceOnePlane(tet, set_VF, correct_normal)));
    auto calculated_volume = getAnalyticVolume(tet, correct_reconstruction[0]);
    EXPECT_NEAR(calculated_volume / tet_volume, set_VF, 1.0e-14)
        << correct_reconstruction << " " << set_VF << std::endl;
  }
}

}  // namespace
