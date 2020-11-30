// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/interface_reconstruction_methods/volume_fraction_matching.h"

#include <algorithm>
#include <cmath>
#include <float.h>
#include <random>

#include "gtest/gtest.h"

#include "irl/generic_cutting/recursive_simplex_cutting/lookup_tables.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/geometry/polyhedrons/hexahedron.h"
#include "irl/helpers/helper.h"
#include "irl/helpers/mymath.h"
#include "irl/planar_reconstruction/planar_localizer.h"
#include "irl/planar_reconstruction/planar_reconstruction.h"
#include "irl/planar_reconstruction/planar_separator.h"

namespace {

using namespace IRL;

TEST(VolumeFractionMatching,
     setDistanceToMatchVolumeFractionPartialFill_RectangularCuboid) {
  double volume_fraction = global_constants::VF_HIGH + DBL_EPSILON;
  PlanarSeparator reconstruction_outer =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 1.0, 0.0), 0.2));
  setDistanceToMatchVolumeFraction(unit_cell, volume_fraction,
                                   &reconstruction_outer);
  EXPECT_DOUBLE_EQ(reconstruction_outer[0].distance(),
                   global_constants::ARBITRARILY_LARGE_DISTANCE);

  volume_fraction = global_constants::VF_LOW - DBL_EPSILON;
  setDistanceToMatchVolumeFraction(unit_cell, volume_fraction,
                                   &reconstruction_outer);
  EXPECT_DOUBLE_EQ(reconstruction_outer[0].distance(),
                   -global_constants::ARBITRARILY_LARGE_DISTANCE);

  std::random_device
      rd; // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd()); // Use the 64-bit Mersenne Twister 19937
                             // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_distance(-0.5, 0.5);
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    double nx = random_normal(eng);
    double ny = random_normal(eng);
    double nz = random_normal(eng);
    double dist = random_distance(eng);
    Plane plane(Normal::normalized(nx, ny, nz), dist);
    PlanarSeparator reconstruction = PlanarSeparator::fromOnePlane(plane);
    double cut_volume = getVolumeFraction(unit_cell, reconstruction);
    setDistanceToMatchVolumeFractionPartialFill(unit_cell, cut_volume,
                                                &reconstruction);
    double new_cut_volume = getVolumeFraction(unit_cell, reconstruction);
    EXPECT_NEAR(new_cut_volume, cut_volume, 1.0e-14);
  }
  std::uniform_real_distribution<double> random_distance_separation(1.0e-6,
                                                                    0.20);
  std::uniform_real_distribution<double> random_distance_shift(-0.15, 0.15);
  std::uniform_real_distribution<double> random_tolerance(1.0e-14, 1.0e-10);
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    double nx = random_normal(eng);
    double ny = random_normal(eng);
    double nz = random_normal(eng);
    double dist = random_distance(eng);
    double dist_separation = random_distance_separation(eng);
    double dist_shift = random_distance_shift(eng);

    double tolerance = random_tolerance(eng);
    Plane plane_0(Normal::normalized(nx, ny, nz), dist);
    Plane plane_1(Normal::normalized(-nx, -ny, -nz),
                  dist + std::copysign(dist_separation, dist));
    PlanarSeparator reconstruction =
        PlanarSeparator::fromTwoPlanes(plane_0, plane_1, 1.0);
    if (dist < 0.0) {
      reconstruction.flipCutting();
    }
    double cut_volume = getVolumeFraction(unit_cell, reconstruction);
    SmallVector<double, 2> distances;
    distances.resize(2);
    distances[0] = plane_0.distance() + dist_shift;
    distances[1] = plane_1.distance() - dist_shift;
    reconstruction.setDistances(distances);
    setDistanceToMatchVolumeFractionPartialFill(unit_cell, cut_volume,
                                                &reconstruction, tolerance);
    double new_cut_volume = getVolumeFraction(unit_cell, reconstruction);
    EXPECT_NEAR(new_cut_volume, cut_volume, tolerance);
  }
  PlanarSeparator all_liquid = PlanarSeparator::fromTwoPlanes(
      Plane(Normal(0.0, -1.0, 0.0), 1000.0),
      Plane(Normal(0.0, 1.0, 0.0), 10000.0), 1.0);
  setDistanceToMatchVolumeFraction(unit_cell, 1.0, &all_liquid);
  EXPECT_EQ(all_liquid.getNumberOfPlanes(), 1);
}

TEST(VolumeFractionMatching,
     setDistanceToMatchVolumeFractionPartialFill_Hexahedron) {
  auto hex_cell =
      Hexahedron({Pt(1.0, 0.0, 0.0), Pt(0.75, 1.0, 0.0), Pt(0.75, 1.0, 1.0),
                  Pt(1.0, 0.0, 1.0), Pt(0.0, 0.0, 0.0), Pt(0.25, 1.0, 0.0),
                  Pt(0.25, 1.0, 1.0), Pt(0.0, 0.0, 1.0)});

  const auto hex_centroid = hex_cell.calculateCentroid();

  double volume_fraction = global_constants::VF_HIGH + DBL_EPSILON;
  PlanarSeparator reconstruction_outer =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 1.0, 0.0), 0.2));
  setDistanceToMatchVolumeFraction(hex_cell, volume_fraction,
                                   &reconstruction_outer, 1.0e-15);
  EXPECT_DOUBLE_EQ(reconstruction_outer[0].distance(),
                   global_constants::ARBITRARILY_LARGE_DISTANCE);

  volume_fraction = global_constants::VF_LOW - DBL_EPSILON;
  setDistanceToMatchVolumeFraction(hex_cell, volume_fraction,
                                   &reconstruction_outer, 1.0e-15);
  EXPECT_DOUBLE_EQ(reconstruction_outer[0].distance(),
                   -global_constants::ARBITRARILY_LARGE_DISTANCE);

  std::random_device
      rd; // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd()); // Use the 64-bit Mersenne Twister 19937
                             // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(global_constants::VF_LOW,
                                                   global_constants::VF_HIGH);
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    double nx = random_normal(eng);
    double ny = random_normal(eng);
    double nz = random_normal(eng);
    const auto normal = Normal::Normal::normalized(nx, ny, nz);
    Plane plane(normal, normal * hex_centroid);
    PlanarSeparator reconstruction = PlanarSeparator::fromOnePlane(plane);
    double cut_volume = getVolumeFraction(hex_cell, reconstruction);
    setDistanceToMatchVolumeFractionPartialFill(hex_cell, cut_volume,
                                                &reconstruction, 1.0e-15);
    double new_cut_volume = getVolumeFraction(hex_cell, reconstruction);
    EXPECT_NEAR(new_cut_volume, cut_volume, 1.0e-14);
  }
}

TEST(VolumeFractionMatching, setDistanceToMatchVolumeFractionPartialFill_Tet) {
  auto tet_cell = Tet({Pt(1.0, 0.0, -0.5), Pt(1.0, 0.5, 0.0), Pt(1.0, 0.0, 0.5),
                       Pt(0.0, 0.0, 0.0)});

  double vf = global_constants::VF_HIGH + DBL_EPSILON;
  PlanarSeparator reconstruction_outer =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 1.0, 0.0), 0.2));
  setDistanceToMatchVolumeFraction(tet_cell, vf, &reconstruction_outer,
                                   1.0e-15);
  EXPECT_DOUBLE_EQ(reconstruction_outer[0].distance(),
                   global_constants::ARBITRARILY_LARGE_DISTANCE);

  vf = global_constants::VF_LOW - DBL_EPSILON;
  setDistanceToMatchVolumeFraction(tet_cell, vf, &reconstruction_outer,
                                   1.0e-15);
  EXPECT_DOUBLE_EQ(reconstruction_outer[0].distance(),
                   -global_constants::ARBITRARILY_LARGE_DISTANCE);

  std::random_device
      rd; // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd()); // Use the 64-bit Mersenne Twister 19937
                             // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> volume_fraction(0.0, 1.0);
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    double nx = random_normal(eng);
    double ny = random_normal(eng);
    double nz = random_normal(eng);
    double dist = 0.0;
    vf = volume_fraction(eng);
    Plane plane(Normal::normalized(nx, ny, nz), dist);
    PlanarSeparator reconstruction = PlanarSeparator::fromOnePlane(plane);
    setDistanceToMatchVolumeFractionPartialFill(tet_cell, vf, &reconstruction,
                                                1.0e-15);
    double new_cut_volume = getVolumeFraction(tet_cell, reconstruction);
    EXPECT_NEAR(new_cut_volume, vf, 1.0e-14);
  }
}

TEST(VolumeFractionMatching, setGroupDistanceToMatchVolumeFraction) {
  // Four Vertically stratified materials
  PlanarSeparator separators[4];
  PlanarSeparatorPathGroup path_group;
  for (UnsignedIndex_t n = 0; n < 4; ++n) {
    separators[n] =
        PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 1.0, 0.0), 0.0));
    PlanarSeparatorPath tmp_path(&separators[n]);
    path_group.addPlanarSeparatorPath(tmp_path, n);
  }
  path_group.setPriorityOrder({0, 1, 2, 3});
  std::array<double, 4> volume_fraction{{0.25, 0.25, 0.25, 0.25}};
  setGroupDistanceToMatchVolumeFraction(unit_cell, volume_fraction, &path_group,
                                        1.0e-15);
  EXPECT_NEAR(separators[0][0].distance(), -0.25, 1.0e-15);
  EXPECT_NEAR(separators[1][0].distance(), 0.0, 1.0e-15);
  EXPECT_NEAR(separators[2][0].distance(), 0.25, 1.0e-15);
  EXPECT_NEAR(separators[3][0].normal()[0], 0.0, 1.0e-15);
  EXPECT_NEAR(separators[3][0].normal()[1], 0.0, 1.0e-15);
  EXPECT_NEAR(separators[3][0].normal()[2], 0.0, 1.0e-15);
  EXPECT_NEAR(separators[3][0].distance(), 1.0, 1.0e-15);

  // Ability to correctly set single phase
  path_group.setPriorityOrder({0});
  std::array<double, 1> single_phase{{1.0}};
  setGroupDistanceToMatchVolumeFraction(unit_cell, single_phase, &path_group,
                                        1.0e-15);
  EXPECT_NEAR(separators[0][0].normal()[0], 0.0, 1.0e-15);
  EXPECT_NEAR(separators[0][0].normal()[1], 0.0, 1.0e-15);
  EXPECT_NEAR(separators[0][0].normal()[2], 0.0, 1.0e-15);
  EXPECT_NEAR(separators[0][0].distance(), 1.0, 1.0e-15);

  // TODO: more complicated arrangment of 4 materials.
}

} // namespace
