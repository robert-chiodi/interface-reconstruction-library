// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/moments/volume_moments.h"

#include <random>

#include "gtest/gtest.h"

#include "src/geometry/general/rotations.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/helpers/geometric_cutting_helpers.h"
#include "src/helpers/helper.h"
#include "src/moments/separated_volume_moments.h"
#include "src/moments/volume_moments_and_normal.h"

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

}  // namespace
