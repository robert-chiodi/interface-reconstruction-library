// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/interface_reconstruction_methods/advected_plane_reconstruction.h"
#include "src/planar_reconstruction/planar_separator.h"

#include <random>

#include "gtest/gtest.h"
#include "src/generic_cutting/generic_cutting.h"

#include "src/helpers/geometric_cutting_helpers.h"
#include "src/moments/volume_moments.h"

namespace {

using namespace IRL;

TEST(AdvectedPlaneReconstruction, AdvectedPlaneReconstruction) {
  // Construct list that has been previously tested in
  // TEST(Partitioning, PartitionByNormalKMeansHarder)
  ListedVolumeMoments<VolumeMomentsAndNormal> list;
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(-0.2, 0.18, 0.0)),
                                 Normal(-0.2, 1.0, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(0.2, -0.18, 0.0)),
                                 Normal(0.2, -1.0, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(0.2, 0.18, 0.0)),
                                 Normal(0.2, 1.0, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(-0.2, -0.18, 0.0)),
                                 Normal(-0.2, -1.0, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.1, Pt(0.0, -0.25, 0.0)),
                                 Normal(0.0, -1.0, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.1, Pt(0.0, 0.25, 0.0)),
                                 Normal(0.0, 1.0, 0.0));

  // Normalize all normals.
  for (auto& element : list) {
    element.normal().normalize();
  }
  list.multiplyByVolume();

  // Set expected PlanarSeparator to create neighborhood moments from.
  PlanarSeparator correct_separator = PlanarSeparator::fromTwoPlanes(
      Plane(Normal(0.0, 1.0, 0.0), -0.25), Plane(Normal(0.0, -1.0, 0.0), -0.25),
      -1.0);
  R2PNeighborhood<RectangularCuboid> neighborhood;
  std::array<RectangularCuboid, 27> cells;
  std::array<SeparatedMoments<VolumeMoments>, 27> moments;
  for (UnsignedIndex_t i = 0; i < 3; ++i) {
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      for (UnsignedIndex_t k = 0; k < 3; ++k) {
        cells[i * 9 + j * 3 + k] = unit_cell;
        cells[i * 9 + j * 3 + k].shift(static_cast<double>(i) - 1.0,
                                       static_cast<double>(j) - 1.0,
                                       static_cast<double>(k) - 1.0);
        moments[i * 9 + j * 3 + k] =
            getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
                cells[i * 9 + j * 3 + k], correct_separator);
        neighborhood.addMember(&cells[i * 9 + j * 3 + k],
                               &moments[i * 9 + j * 3 + k]);
      }
    }
  }
  neighborhood.setCenterOfStencil(13);

  auto computed_separator =
      AdvectedPlaneReconstruction::solve(list, neighborhood);

  ASSERT_EQ(computed_separator.getNumberOfPlanes(), 2);
  EXPECT_NEAR(computed_separator[0].normal()[0], 0.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[0].normal()[1], 1.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[0].normal()[2], 0.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[0].distance(), -0.25, 1.0e-8);
  EXPECT_NEAR(computed_separator[1].normal()[0], 0.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[1].normal()[1], -1.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[1].normal()[2], 0.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[1].distance(), -0.25, 1.0e-8);
  EXPECT_EQ(computed_separator.isFlipped(), true);
}

TEST(AdvectedPlaneReconstruction, AdvectedPlaneReconstructionShifted) {
  // Construct list that has been previously tested in
  // TEST(Partitioning, PartitionByNormalKMeansHarder)
  ListedVolumeMoments<VolumeMomentsAndNormal> list;
  list += VolumeMomentsAndNormal(
      VolumeMoments(10.0 * 0.5, 10.0 * Pt(-0.2, 0.18, 0.0)),
      Normal(-0.2, 1.0, 0.0));
  list += VolumeMomentsAndNormal(
      VolumeMoments(10.0 * 0.5, 10.0 * Pt(0.2, -0.18, 0.0)),
      Normal(0.2, -1.0, 0.0));
  list += VolumeMomentsAndNormal(
      VolumeMoments(10.0 * 0.5, 10.0 * Pt(0.2, 0.18, 0.0)),
      Normal(0.2, 1.0, 0.0));
  list += VolumeMomentsAndNormal(
      VolumeMoments(10.0 * 0.5, 10.0 * Pt(-0.2, -0.18, 0.0)),
      Normal(-0.2, -1.0, 0.0));
  list += VolumeMomentsAndNormal(
      VolumeMoments(10.0 * 0.1, 10.0 * Pt(0.0, -0.25, 0.0)),
      Normal(0.0, -1.0, 0.0));
  list += VolumeMomentsAndNormal(
      VolumeMoments(10.0 * 0.1, 10.0 * Pt(0.0, 0.25, 0.0)),
      Normal(0.0, 1.0, 0.0));

  // Normalize all normals and shift centroids.
  for (auto& element : list) {
    element.volumeMoments().centroid() += Pt(-100.0, -100.0, -100.0);
    element.normal().normalize();
  }
  list.multiplyByVolume();

  // Set expected PlanarSeparator to create neighborhood moments from.
  PlanarSeparator correct_separator = PlanarSeparator::fromTwoPlanes(
      Plane(Normal(0.0, 1.0, 0.0), -100.25),
      Plane(Normal(0.0, -1.0, 0.0), 99.75), -1.0);
  RectangularCuboid middle_cell = RectangularCuboid::fromBoundingPts(
      Pt(-105.0, -105.0, -105.0), Pt(-95.0, -95.0, -95.0));
  R2PNeighborhood<RectangularCuboid> neighborhood;
  std::array<RectangularCuboid, 27> cells;
  std::array<SeparatedMoments<VolumeMoments>, 27> moments;
  for (UnsignedIndex_t i = 0; i < 3; ++i) {
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      for (UnsignedIndex_t k = 0; k < 3; ++k) {
        cells[i * 9 + j * 3 + k] = middle_cell;
        cells[i * 9 + j * 3 + k].shift(10.0 * (static_cast<double>(i) - 1.0),
                                       10.0 * (static_cast<double>(j) - 1.0),
                                       10.0 * (static_cast<double>(k) - 1.0));
        moments[i * 9 + j * 3 + k] =
            getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
                cells[i * 9 + j * 3 + k], correct_separator);
        neighborhood.addMember(&cells[i * 9 + j * 3 + k],
                               &moments[i * 9 + j * 3 + k]);
      }
    }
  }
  neighborhood.setCenterOfStencil(13);
  auto computed_separator =
      AdvectedPlaneReconstruction::solve(list, neighborhood);
  ASSERT_EQ(computed_separator.getNumberOfPlanes(), 2);
  EXPECT_NEAR(computed_separator[0].normal()[0], 0.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[0].normal()[1], 1.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[0].normal()[2], 0.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[0].distance(), -100.25, 1.0e-8);
  EXPECT_NEAR(computed_separator[1].normal()[0], 0.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[1].normal()[1], -1.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[1].normal()[2], 0.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[1].distance(), 99.75, 1.0e-8);
  EXPECT_EQ(computed_separator.isFlipped(), true);
}

}  // namespace
