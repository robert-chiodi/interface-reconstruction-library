// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/moments/volume_moments.h"

#include <random>

#include "gtest/gtest.h"

#include "irl/geometry/general/rotations.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/helpers/geometric_cutting_helpers.h"
#include "irl/helpers/helper.h"
#include "irl/moments/separated_volume_moments.h"
#include "irl/moments/volume_moments_and_normal.h"

namespace {

using namespace IRL;

TEST(VolumeMoments, Volume) {
  Volume volume;
  EXPECT_DOUBLE_EQ(volume, 0.0);
  volume = Volume(3.8);
  EXPECT_DOUBLE_EQ(volume, 3.8);
  volume = Volume::calculateMoments(&unit_cell);
  EXPECT_NEAR(volume, 1.0, 1.0e-15);
  volume.multiplyByVolume();
  EXPECT_NEAR(volume, 1.0, 1.0e-15);
  volume.normalizeByVolume();
  EXPECT_NEAR(volume, 1.0, 1.0e-15);
  volume += Volume(2.0);
  EXPECT_NEAR(volume, 3.0, 1.0e-15);
  volume *= 2.0;
  EXPECT_NEAR(volume, 6.0, 1.0e-15);
  volume = -1.0;
  EXPECT_DOUBLE_EQ(volume, -1.0);
}

TEST(VolumeMoments, VolumeMoments) {
  // Setting, Accessing
  VolumeMoments init_moment;
  EXPECT_DOUBLE_EQ(init_moment.volume(), 0.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[0], 0.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[1], 0.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[2], 0.0);

  VolumeMoments nonzero_moment(0.25, Pt(-0.5, 1.0, -8.0));
  EXPECT_DOUBLE_EQ(nonzero_moment.volume(), 0.25);
  EXPECT_DOUBLE_EQ(nonzero_moment.centroid()[0], -0.5);
  EXPECT_DOUBLE_EQ(nonzero_moment.centroid()[1], 1.0);
  EXPECT_DOUBLE_EQ(nonzero_moment.centroid()[2], -8.0);

  // Operators
  init_moment += nonzero_moment;
  EXPECT_DOUBLE_EQ(init_moment.volume(), 0.25);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[0], -0.5);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[1], 1.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[2], -8.0);

  init_moment = init_moment - nonzero_moment;
  EXPECT_DOUBLE_EQ(init_moment.volume(), 0.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[0], 0.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[1], 0.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[2], 0.0);

  init_moment = init_moment + nonzero_moment;
  EXPECT_DOUBLE_EQ(init_moment.volume(), 0.25);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[0], -0.5);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[1], 1.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[2], -8.0);

  // Volume normalization and multiplication
  init_moment.multiplyByVolume();
  EXPECT_DOUBLE_EQ(init_moment.volume(), 0.25);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[0], -0.5 * 0.25);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[1], 1.0 * 0.25);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[2], -8.0 * 0.25);

  init_moment.normalizeByVolume();
  EXPECT_DOUBLE_EQ(init_moment.volume(), 0.25);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[0], -0.5);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[1], 1.0);
  EXPECT_DOUBLE_EQ(init_moment.centroid()[2], -8.0);
}

TEST(VolumeMoments, SeparatedMoments) {
  VolumeMoments liquid_moments(0.25, Pt(-0.5, 1.0, -8.0));
  VolumeMoments gas_moments(0.80, Pt(0.5, -1.0, 8.0));

  // Construction and access
  SeparatedMoments<VolumeMoments> blank_moments;
  SeparatedMoments<VolumeMoments> nonzero_moments(liquid_moments, gas_moments);
  EXPECT_DOUBLE_EQ(nonzero_moments[0].volume(), 0.25);
  EXPECT_DOUBLE_EQ(nonzero_moments[0].centroid()[0], -0.5);
  EXPECT_DOUBLE_EQ(nonzero_moments[0].centroid()[1], 1.0);
  EXPECT_DOUBLE_EQ(nonzero_moments[0].centroid()[2], -8.0);
  EXPECT_DOUBLE_EQ(nonzero_moments[1].volume(), 0.80);
  EXPECT_DOUBLE_EQ(nonzero_moments[1].centroid()[0], 0.5);
  EXPECT_DOUBLE_EQ(nonzero_moments[1].centroid()[1], -1.0);
  EXPECT_DOUBLE_EQ(nonzero_moments[1].centroid()[2], 8.0);

  // Operators
  blank_moments += nonzero_moments;
  EXPECT_DOUBLE_EQ(blank_moments[0].volume(), 0.25);
  EXPECT_DOUBLE_EQ(blank_moments[0].centroid()[0], -0.5);
  EXPECT_DOUBLE_EQ(blank_moments[0].centroid()[1], 1.0);
  EXPECT_DOUBLE_EQ(blank_moments[0].centroid()[2], -8.0);
  EXPECT_DOUBLE_EQ(blank_moments[1].volume(), 0.80);
  EXPECT_DOUBLE_EQ(blank_moments[1].centroid()[0], 0.5);
  EXPECT_DOUBLE_EQ(blank_moments[1].centroid()[1], -1.0);
  EXPECT_DOUBLE_EQ(blank_moments[1].centroid()[2], 8.0);

  // Multiplication and normalization of each phase by their volumes
  nonzero_moments.multiplyByVolume();
  EXPECT_DOUBLE_EQ(nonzero_moments[0].volume(), 0.25);
  EXPECT_DOUBLE_EQ(nonzero_moments[0].centroid()[0], -0.5 * 0.25);
  EXPECT_DOUBLE_EQ(nonzero_moments[0].centroid()[1], 1.0 * 0.25);
  EXPECT_DOUBLE_EQ(nonzero_moments[0].centroid()[2], -8.0 * 0.25);
  EXPECT_DOUBLE_EQ(nonzero_moments[1].volume(), 0.80);
  EXPECT_DOUBLE_EQ(nonzero_moments[1].centroid()[0], 0.5 * 0.80);
  EXPECT_DOUBLE_EQ(nonzero_moments[1].centroid()[1], -1.0 * 0.80);
  EXPECT_DOUBLE_EQ(nonzero_moments[1].centroid()[2], 8.0 * 0.80);

  nonzero_moments.normalizeByVolume();
  EXPECT_DOUBLE_EQ(nonzero_moments[0].volume(), 0.25);
  EXPECT_DOUBLE_EQ(nonzero_moments[0].centroid()[0], -0.5);
  EXPECT_DOUBLE_EQ(nonzero_moments[0].centroid()[1], 1.0);
  EXPECT_DOUBLE_EQ(nonzero_moments[0].centroid()[2], -8.0);
  EXPECT_DOUBLE_EQ(nonzero_moments[1].volume(), 0.80);
  EXPECT_DOUBLE_EQ(nonzero_moments[1].centroid()[0], 0.5);
  EXPECT_DOUBLE_EQ(nonzero_moments[1].centroid()[1], -1.0);
  EXPECT_DOUBLE_EQ(nonzero_moments[1].centroid()[2], 8.0);
}

TEST(VolumeMoments,
     SeparatedVolumeMoments_getSeparatedVolumeMomentsFromComplement) {
  // Testing finding of complement of a cube
  RectangularCuboid unit_cube = RectangularCuboid::fromBoundingPts(
      Pt(-0.5, -0.5, -0.5), Pt(0.5, 0.5, 0.5));
  VolumeMoments given_moments_for_cube(0.3, Pt(0.0, -0.5 + 0.15, 0.0));
  given_moments_for_cube.multiplyByVolume();
  SeparatedMoments<VolumeMoments> full_from_cube =
      SeparatedMoments<VolumeMoments>::fillWithComplementMoments(
          given_moments_for_cube, unit_cube, false);
  full_from_cube.normalizeByVolume();
  EXPECT_DOUBLE_EQ(full_from_cube[0].volume(), 0.3);
  EXPECT_DOUBLE_EQ(full_from_cube[0].centroid()[0], 0.0);
  EXPECT_DOUBLE_EQ(full_from_cube[0].centroid()[1], -0.35);
  EXPECT_DOUBLE_EQ(full_from_cube[0].centroid()[2], 0.0);
  EXPECT_DOUBLE_EQ(full_from_cube[1].volume(), 0.7);
  EXPECT_DOUBLE_EQ(full_from_cube[1].centroid()[0], 0.0);
  EXPECT_DOUBLE_EQ(full_from_cube[1].centroid()[1], 0.15);
  EXPECT_DOUBLE_EQ(full_from_cube[1].centroid()[2], 0.0);

  // Testing finding of complement of a tet
  Tet tet({Pt(-1.0, -1.0, -1.0), Pt(-1.0, 0.0, 0.0), Pt(0.0, 0.0, -1.0),
           Pt(0.0, 0.0, 0.0)});
  double tet_vol = tet.calculateVolume();
  Pt tet_centroid = tet.calculateCentroid();
  VolumeMoments given_moments_for_tet(0.1, Pt(-0.25, -0.125, -0.125));
  given_moments_for_tet.multiplyByVolume();
  SeparatedMoments<VolumeMoments> full_from_tet =
      SeparatedMoments<VolumeMoments>::fillWithComplementMoments(
          given_moments_for_tet, tet, false);
  full_from_tet.normalizeByVolume();
  EXPECT_DOUBLE_EQ(full_from_tet[0].volume(), 0.1);
  EXPECT_DOUBLE_EQ(full_from_tet[0].centroid()[0], -0.25);
  EXPECT_DOUBLE_EQ(full_from_tet[0].centroid()[1], -0.125);
  EXPECT_DOUBLE_EQ(full_from_tet[0].centroid()[2], -0.125);
  EXPECT_DOUBLE_EQ(full_from_tet[1].volume(), 1.0 / 6.0 - 0.1);
  EXPECT_DOUBLE_EQ(full_from_tet[1].centroid()[0],
                   (tet_vol * tet_centroid[0] - 0.1 * -0.25) / (tet_vol - 0.1));
  EXPECT_DOUBLE_EQ(
      full_from_tet[1].centroid()[1],
      (tet_vol * tet_centroid[1] - 0.1 * -0.125) / (tet_vol - 0.1));
  EXPECT_DOUBLE_EQ(
      full_from_tet[1].centroid()[2],
      (tet_vol * tet_centroid[2] - 0.1 * -0.125) / (tet_vol - 0.1));

  // Check extremes (full liquid or gas) using the tet
  given_moments_for_tet.volume() = (0.0);
  given_moments_for_tet.multiplyByVolume();
  full_from_tet = SeparatedMoments<VolumeMoments>::fillWithComplementMoments(
      given_moments_for_tet, tet, false);
  full_from_tet.normalizeByVolume();
  EXPECT_DOUBLE_EQ(full_from_tet[0].volume(), 0.0);
  EXPECT_DOUBLE_EQ(full_from_tet[0].centroid()[0], 0.0);
  EXPECT_DOUBLE_EQ(full_from_tet[0].centroid()[1], 0.0);
  EXPECT_DOUBLE_EQ(full_from_tet[0].centroid()[2], 0.0);
  EXPECT_DOUBLE_EQ(full_from_tet[1].volume(), tet_vol);
  EXPECT_DOUBLE_EQ(full_from_tet[1].centroid()[0], tet_centroid[0]);
  EXPECT_DOUBLE_EQ(full_from_tet[1].centroid()[1], tet_centroid[1]);
  EXPECT_DOUBLE_EQ(full_from_tet[1].centroid()[2], tet_centroid[2]);

  given_moments_for_tet.volume() = (tet_vol);
  given_moments_for_tet.centroid() = (tet_centroid);
  given_moments_for_tet.multiplyByVolume();
  full_from_tet = SeparatedMoments<VolumeMoments>::fillWithComplementMoments(
      given_moments_for_tet, tet, false);
  full_from_tet.normalizeByVolume();
  EXPECT_DOUBLE_EQ(full_from_tet[0].volume(), tet_vol);
  EXPECT_DOUBLE_EQ(full_from_tet[0].centroid()[0], tet_centroid[0]);
  EXPECT_DOUBLE_EQ(full_from_tet[0].centroid()[1], tet_centroid[1]);
  EXPECT_DOUBLE_EQ(full_from_tet[0].centroid()[2], tet_centroid[2]);
  EXPECT_DOUBLE_EQ(full_from_tet[1].volume(), 0.0);
  EXPECT_DOUBLE_EQ(full_from_tet[1].centroid()[0], 0.0);
  EXPECT_DOUBLE_EQ(full_from_tet[1].centroid()[1], 0.0);
  EXPECT_DOUBLE_EQ(full_from_tet[1].centroid()[2], 0.0);
}

TEST(VolumeMoments, VolumeMomentsAndNormal) {
  // Setting, Accessing
  VolumeMomentsAndNormal init_moment;
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().volume(), 0.0);
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().centroid()[0], 0.0);
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().centroid()[1], 0.0);
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().centroid()[2], 0.0);
  EXPECT_DOUBLE_EQ(init_moment.normal()[0], 0.0);
  EXPECT_DOUBLE_EQ(init_moment.normal()[1], 0.0);
  EXPECT_DOUBLE_EQ(init_moment.normal()[2], 0.0);

  VolumeMomentsAndNormal nonzero_moment(
      VolumeMoments(0.25, Pt(-0.5, 1.0, -8.0)), Normal(0.0, -1.0, 0.0));
  EXPECT_DOUBLE_EQ(nonzero_moment.volumeMoments().volume(), 0.25);
  EXPECT_DOUBLE_EQ(nonzero_moment.volumeMoments().centroid()[0], -0.5);
  EXPECT_DOUBLE_EQ(nonzero_moment.volumeMoments().centroid()[1], 1.0);
  EXPECT_DOUBLE_EQ(nonzero_moment.volumeMoments().centroid()[2], -8.0);
  EXPECT_DOUBLE_EQ(nonzero_moment.normal()[0], 0.0);
  EXPECT_DOUBLE_EQ(nonzero_moment.normal()[1], -1.0);
  EXPECT_DOUBLE_EQ(nonzero_moment.normal()[2], 0.0);

  // Operators
  init_moment += nonzero_moment;
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().volume(), 0.25);
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().centroid()[0], -0.5);
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().centroid()[1], 1.0);
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().centroid()[2], -8.0);
  EXPECT_DOUBLE_EQ(init_moment.normal()[0], 0.0);
  EXPECT_DOUBLE_EQ(init_moment.normal()[1], -1.0);
  EXPECT_DOUBLE_EQ(init_moment.normal()[2], 0.0);

  // Volume normalization and multiplication
  init_moment.multiplyByVolume();
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().volume(), 0.25);
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().centroid()[0], -0.5 * 0.25);
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().centroid()[1], 1.0 * 0.25);
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().centroid()[2], -8.0 * 0.25);
  EXPECT_DOUBLE_EQ(init_moment.normal()[0], 0.0);
  EXPECT_DOUBLE_EQ(init_moment.normal()[1], -0.25);
  EXPECT_DOUBLE_EQ(init_moment.normal()[2], 0.0);

  init_moment.normalizeByVolume();
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().volume(), 0.25);
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().centroid()[0], -0.5);
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().centroid()[1], 1.0);
  EXPECT_DOUBLE_EQ(init_moment.volumeMoments().centroid()[2], -8.0);
  EXPECT_DOUBLE_EQ(init_moment.normal()[0], 0.0);
  EXPECT_DOUBLE_EQ(init_moment.normal()[1], -1.0);
  EXPECT_DOUBLE_EQ(init_moment.normal()[2], 0.0);
}

TEST(VolumeMoments, GeneralMoments3D) {
  auto cube =
      RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(1.0, 1.0, 1.0));
  auto mom = cube.calculateGeneralMoments<2>();
  std::array<double, 10> correct{
      {1.0, 0.5, 0.5, 0.5, 1.0 / 3.0, 0.25, 0.25, 1.0 / 3.0, 0.25, 1.0 / 3.0}};
  for (UnsignedIndex_t i = 0; i < mom.size(); ++i) {
    EXPECT_DOUBLE_EQ(mom[i], correct[i]);
  }

  auto tet = Tet({
      Pt(2.0, 0.0, 0.0),
      Pt(10.0, 3.0, 2.0),
      Pt(10.0, 5.0, 0.0),
      Pt(10.0, 4.0, -3.0),
  });
  // Correct generated using mathematica to do the integrations
  correct = std::array<double, 10>{{10.6666666666666, 85.3333333333333, 32.0,
                                    -2.666666666666667, 708.2666666666666,
                                    268.8, -22.4, 103.4666666666667, -9.6,
                                    7.466666666666666}};
  mom = tet.calculateGeneralMoments<2>();
  for (UnsignedIndex_t i = 0; i < mom.size(); ++i) {
    EXPECT_NEAR(mom[i], correct[i], 1.0e-12);
  }
}

}  // namespace
