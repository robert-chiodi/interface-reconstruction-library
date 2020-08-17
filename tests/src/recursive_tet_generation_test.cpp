// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/generic_cutting/generic_cutting.h"

#include "gtest/gtest.h"

#include <array>

#include "src/geometry/polyhedrons/tet.h"
#include "src/helpers/geometric_cutting_helpers.h"

namespace {

using namespace IRL;

TEST(RecursiveTetGeneration, TetKeepsPositiveSign) {
  Tet positive_tet({Pt(1.0, 0.0, -0.5), Pt(1.0, 1.0, 0.0), Pt(1.0, 0.0, 0.5),
                    Pt(0.0, 0.0, 0.0)});
  auto tet_volume = positive_tet.calculateVolume();
  ASSERT_TRUE(tet_volume > 0.0);

  std::array<double, 4> distance_to_vertices;

  // Now test all the cutting cases and make sure resulting volume for
  // under/over plane is positive.
  PlanarSeparator planar_reconstruction =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 0.0), 0.0));

  // case 0
  Normal plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 2.0);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 0);
  auto separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                       RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 1
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 1);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 2
  plane_normal = Normal(0.0, 1.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 2);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 3
  plane_normal = Normal::normalized(1.0, 0.0, -1.0);
  planar_reconstruction[0] =
      Plane(plane_normal, positive_tet[2] * plane_normal * 1.1);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 3);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 4
  plane_normal = Normal::normalized(0.0, 0.0, 1.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 4);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 5
  plane_normal = Normal::normalized(1.0, -1.0, 0.0);
  planar_reconstruction[0] =
      Plane(plane_normal, positive_tet[2] * plane_normal * 0.9);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 5);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 6
  plane_normal = Normal::normalized(1.0, 0.0, 1.0);
  planar_reconstruction[0] =
      Plane(plane_normal, positive_tet[1] * plane_normal * 0.9);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 6);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 7
  plane_normal = Normal::normalized(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.5);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 7);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 15
  plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -2.0);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 15);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 14
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 14);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 13
  plane_normal = Normal(0.0, 1.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 13);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 12
  plane_normal = Normal::normalized(1.0, 0.0, -1.0);
  planar_reconstruction[0] =
      Plane(-plane_normal, -positive_tet[2] * plane_normal * 1.1);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 12);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 11
  plane_normal = Normal::normalized(0.0, 0.0, 1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 11);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 10
  plane_normal = Normal::normalized(1.0, -1.0, 0.0);
  planar_reconstruction[0] =
      Plane(-plane_normal, -positive_tet[2] * plane_normal * 0.9);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 10);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 9
  plane_normal = Normal::normalized(1.0, 0.0, 1.0);
  planar_reconstruction[0] =
      Plane(-plane_normal, -positive_tet[1] * plane_normal * 0.9);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 9);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 8
  plane_normal = Normal::normalized(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.5);
  signedDistanceToVertices(positive_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 8);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);
}

TEST(RecursiveTetGeneration, TetKeepsNegativeSign) {
  Tet negative_tet({Pt(1.0, 0.0, -0.5), Pt(1.0, 0.0, 0.5), Pt(1.0, 1.0, 0.0),
                    Pt(0.0, 0.0, 0.0)});
  auto tet_volume = negative_tet.calculateVolume();
  ASSERT_TRUE(tet_volume < 0.0);

  std::array<double, 4> distance_to_vertices;

  // Now test all the cutting cases and make sure resulting volume for
  // under/over plane is positive.
  PlanarSeparator planar_reconstruction =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 0.0), 0.0));

  // case 0
  Normal plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 2.0);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 0);
  auto separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                       RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume(), 1.0e-14);

  // case 1
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 1);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 2
  plane_normal = Normal::normalized(0.0, 0.0, 1.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 2);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 3
  plane_normal = Normal::normalized(1.0, -1.0, 0.0);
  planar_reconstruction[0] =
      Plane(plane_normal, negative_tet[1] * plane_normal * 0.9);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 3);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 4
  plane_normal = Normal(0.0, 1.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 4);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 5
  plane_normal = Normal::normalized(1.0, 0.0, -1.0);
  planar_reconstruction[0] =
      Plane(plane_normal, negative_tet[1] * plane_normal * 1.2);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 5);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 6
  plane_normal = Normal::normalized(1.0, 0.0, 1.0);
  planar_reconstruction[0] =
      Plane(plane_normal, negative_tet[0] * plane_normal * 1.1);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 6);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 7
  plane_normal = Normal::normalized(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.5);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 7);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 15
  plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -2.0);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 15);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[1].volume(), 1.0e-14);

  // case 14
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 14);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 13
  plane_normal = Normal::normalized(0.0, 0.0, 1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 13);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 12
  plane_normal = Normal::normalized(1.0, -1.0, 0.0);
  planar_reconstruction[0] =
      Plane(-plane_normal, -negative_tet[1] * plane_normal * 0.9);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 12);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 11
  plane_normal = Normal(0.0, 1.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 11);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 10
  plane_normal = Normal::normalized(1.0, 0.0, -1.0);
  planar_reconstruction[0] =
      Plane(-plane_normal, -negative_tet[1] * plane_normal * 1.2);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 10);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 9
  plane_normal = Normal::normalized(1.0, 0.0, 1.0);
  planar_reconstruction[0] =
      Plane(-plane_normal, -negative_tet[0] * plane_normal * 1.1);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 9);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 8
  plane_normal = Normal::normalized(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.5);
  signedDistanceToVertices(negative_tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 8);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tet, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tet_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);
}

}  // namespace
