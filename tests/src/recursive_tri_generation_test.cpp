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

#include "src/geometry/polygons/tri.h"
#include "src/helpers/geometric_cutting_helpers.h"

namespace {

using namespace IRL;

TEST(RecursiveTriGeneration, TriKeepsPositiveSign) {
  Tri positive_tri({Pt(1.0, 0.0, -0.5), Pt(1.0, 1.0, 0.0), Pt(1.0, 0.0, 0.5)});
  positive_tri.calculateAndSetPlaneOfExistence();
  auto plane_of_existence = positive_tri.getPlaneOfExistence();
  ASSERT_NEAR(plane_of_existence.normal()[0], 1.0, 1.0e-15);
  ASSERT_NEAR(plane_of_existence.normal()[1], 0.0, 1.0e-15);
  ASSERT_NEAR(plane_of_existence.normal()[2], 0.0, 1.0e-15);
  ASSERT_NEAR(plane_of_existence.distance(), 1.0, 1.0e-15);
  auto tri_volume = positive_tri.calculateVolume();
  ASSERT_TRUE(tri_volume > 0.0);

  std::array<double, 3> distance_to_vertices;

  // Now test all the cutting cases and make sure resulting volume for
  // under/over plane is positive.
  PlanarSeparator planar_reconstruction =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 0.0), 0.0));

  // case 0
  Normal plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 2.0);
  signedDistanceToVertices(positive_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 0);
  auto separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                       RecursiveSimplexCutting>(
      positive_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 1
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(positive_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 1);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 2
  plane_normal = Normal(0.0, 1.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(positive_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 2);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 3
  plane_normal = Normal::normalized(1.0, 0.0, -1.0);
  planar_reconstruction[0] =
      Plane(plane_normal, positive_tri[2] * plane_normal * 1.1);
  signedDistanceToVertices(positive_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 3);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 7
  plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -2.0);
  signedDistanceToVertices(positive_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 7);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 6
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(positive_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 6);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 5
  plane_normal = Normal(0.0, 1.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(positive_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 5);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 4
  plane_normal = Normal::normalized(1.0, 0.0, -1.0);
  planar_reconstruction[0] =
      Plane(-plane_normal, -positive_tri[2] * plane_normal * 1.1);
  signedDistanceToVertices(positive_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 4);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      positive_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);
}

TEST(RecursiveTriGeneration, FlippedTriKeepsPositiveSign) {
  Tri flipped_tri({Pt(1.0, 0.0, -0.5), Pt(1.0, 0.0, 0.5), Pt(1.0, 1.0, 0.0)});
  flipped_tri.calculateAndSetPlaneOfExistence();
  auto plane_of_existence = flipped_tri.getPlaneOfExistence();
  ASSERT_NEAR(plane_of_existence.normal()[0], -1.0, 1.0e-15);
  ASSERT_NEAR(plane_of_existence.normal()[1], 0.0, 1.0e-15);
  ASSERT_NEAR(plane_of_existence.normal()[2], 0.0, 1.0e-15);
  ASSERT_NEAR(plane_of_existence.distance(), -1.0, 1.0e-15);
  auto tri_volume = flipped_tri.calculateVolume();
  ASSERT_TRUE(tri_volume > 0.0);

  std::array<double, 3> distance_to_vertices;

  // Now test all the cutting cases and make sure resulting volume for
  // under/over plane is positive.
  PlanarSeparator planar_reconstruction =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 0.0), 0.0));

  // case 0
  Normal plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 2.0);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 0);
  auto separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                       RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 1
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 1);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 2
  plane_normal = Normal::normalized(0.0, 0.0, 1.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 2);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 3
  plane_normal = Normal::normalized(1.0, -1.0, 0.0);
  planar_reconstruction[0] =
      Plane(plane_normal, flipped_tri[1] * plane_normal * 0.9);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 3);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 7
  plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -2.0);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 7);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 6
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 6);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 5
  plane_normal = Normal::normalized(0.0, 0.0, 1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 5);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 4
  plane_normal = Normal::normalized(1.0, -1.0, 0.0);
  planar_reconstruction[0] =
      Plane(-plane_normal, -flipped_tri[1] * plane_normal * 0.9);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 4);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() > 0.0);
  ASSERT_TRUE(separator_vm[1].volume() > 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);
}

TEST(RecursiveTriGeneration, TriKeepsNegativeSign) {
  Tri negative_tri({Pt(1.0, 0.0, -0.5), Pt(1.0, 1.0, 0.0), Pt(1.0, 0.0, 0.5)});
  negative_tri.calculateAndSetPlaneOfExistence();
  auto plane_of_existence = negative_tri.getPlaneOfExistence();
  negative_tri.setPlaneOfExistence(plane_of_existence.generateFlippedPlane());
  plane_of_existence = negative_tri.getPlaneOfExistence();
  ASSERT_NEAR(plane_of_existence.normal()[0], -1.0, 1.0e-15);
  ASSERT_NEAR(plane_of_existence.normal()[1], 0.0, 1.0e-15);
  ASSERT_NEAR(plane_of_existence.normal()[2], 0.0, 1.0e-15);
  ASSERT_NEAR(plane_of_existence.distance(), -1.0, 1.0e-15);
  auto tri_volume = negative_tri.calculateVolume();
  ASSERT_TRUE(tri_volume < 0.0);

  std::array<double, 3> distance_to_vertices;

  // Now test all the cutting cases and make sure resulting volume for
  // under/over plane is positive.
  PlanarSeparator planar_reconstruction =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 0.0), 0.0));

  // case 0
  Normal plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 2.0);
  signedDistanceToVertices(negative_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 0);
  auto separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                       RecursiveSimplexCutting>(
      negative_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 1
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(negative_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 1);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 2
  plane_normal = Normal(0.0, 1.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(negative_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 2);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 3
  plane_normal = Normal::normalized(1.0, 0.0, -1.0);
  planar_reconstruction[0] =
      Plane(plane_normal, negative_tri[2] * plane_normal * 1.1);
  signedDistanceToVertices(negative_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 3);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 7
  plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -2.0);
  signedDistanceToVertices(negative_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 7);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 6
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(negative_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 6);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 5
  plane_normal = Normal(0.0, 1.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(negative_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 5);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 4
  plane_normal = Normal::normalized(1.0, 0.0, -1.0);
  planar_reconstruction[0] =
      Plane(-plane_normal, -negative_tri[2] * plane_normal * 1.1);
  signedDistanceToVertices(negative_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 4);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      negative_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);
}

TEST(RecursiveTriGeneration, FlippedTriKeepsNegativeSign) {
  Tri flipped_tri({Pt(1.0, 0.0, -0.5), Pt(1.0, 0.0, 0.5), Pt(1.0, 1.0, 0.0)});
  flipped_tri.calculateAndSetPlaneOfExistence();
  auto plane_of_existence = flipped_tri.getPlaneOfExistence();
  flipped_tri.setPlaneOfExistence(plane_of_existence.generateFlippedPlane());
  plane_of_existence = flipped_tri.getPlaneOfExistence();
  ASSERT_NEAR(plane_of_existence.normal()[0], 1.0, 1.0e-15);
  ASSERT_NEAR(plane_of_existence.normal()[1], 0.0, 1.0e-15);
  ASSERT_NEAR(plane_of_existence.normal()[2], 0.0, 1.0e-15);
  ASSERT_NEAR(plane_of_existence.distance(), 1.0, 1.0e-15);
  auto tri_volume = flipped_tri.calculateVolume();
  ASSERT_TRUE(tri_volume < 0.0);

  std::array<double, 3> distance_to_vertices;

  // Now test all the cutting cases and make sure resulting volume for
  // under/over plane is positive.
  PlanarSeparator planar_reconstruction =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 0.0), 0.0));

  // case 0
  Normal plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 2.0);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 0);
  auto separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                       RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 1
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 1);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 2
  plane_normal = Normal::normalized(0.0, 0.0, 1.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 2);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 3
  plane_normal = Normal::normalized(1.0, -1.0, 0.0);
  planar_reconstruction[0] =
      Plane(plane_normal, flipped_tri[1] * plane_normal * 0.9);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 3);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 7
  plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -2.0);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 7);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 6
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 6);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 5
  plane_normal = Normal::normalized(0.0, 0.0, 1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 5);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);

  // case 4
  plane_normal = Normal::normalized(1.0, -1.0, 0.0);
  planar_reconstruction[0] =
      Plane(-plane_normal, -flipped_tri[1] * plane_normal * 0.9);
  signedDistanceToVertices(flipped_tri, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 4);
  separator_vm = getVolumeMoments<SeparatedMoments<VolumeMoments>,
                                  RecursiveSimplexCutting>(
      flipped_tri, planar_reconstruction);
  ASSERT_TRUE(separator_vm[0].volume() < 0.0);
  ASSERT_TRUE(separator_vm[1].volume() < 0.0);
  ASSERT_NEAR(tri_volume, separator_vm[0].volume() + separator_vm[1].volume(),
              1.0e-14);
}

}  // namespace
