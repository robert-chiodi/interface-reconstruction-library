// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/distributions/partition_by_normal_vector.h"

#include <random>

#include "gtest/gtest.h"

#include "src/distributions/k_means.h"
#include "src/geometry/general/normal.h"
#include "src/moments/listed_volume_moments.h"
#include "src/moments/volume_moments_and_normal.h"

namespace {

using namespace IRL;

TEST(Partitioning, PartitionByNormalKMeans) {
  ListedVolumeMoments<VolumeMomentsAndNormal> list;
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(-0.5, -1.0, -2.0)),
                                 Normal::normalized(-1.0, -0.3, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(0.5, 1.0, 2.0)),
                                 Normal::normalized(-1.0, 0.3, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(-3.0, -6.0, -4.0)),
                                 Normal::normalized(1.0, -0.3, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(3.0, 6.0, 4.0)),
                                 Normal::normalized(1.0, 0.3, 0.0));

  list.multiplyByVolume();

  PartitionByNormal<ListedVolumeMoments<VolumeMomentsAndNormal>> partitioned(
      &list);
  partitioned.setup();
  auto iter = KMeans::partition(&partitioned);
  auto moments = partitioned.getPartitionedObjects();
  moments[0].normalizeByVolume();
  moments[1].normalizeByVolume();
  moments[0].normal().normalize();
  moments[1].normal().normalize();
  EXPECT_DOUBLE_EQ(moments[0].volumeMoments().volume(), 1.0);
  EXPECT_DOUBLE_EQ(moments[0].volumeMoments().centroid()[0], 0.0);
  EXPECT_DOUBLE_EQ(moments[0].volumeMoments().centroid()[1], 0.0);
  EXPECT_DOUBLE_EQ(moments[0].volumeMoments().centroid()[2], 0.0);
  EXPECT_DOUBLE_EQ(moments[0].normal()[0], -1.0);
  EXPECT_DOUBLE_EQ(moments[0].normal()[1], 0.0);
  EXPECT_DOUBLE_EQ(moments[0].normal()[2], 0.0);
  EXPECT_DOUBLE_EQ(moments[1].volumeMoments().volume(), 1.0);
  EXPECT_DOUBLE_EQ(moments[1].volumeMoments().centroid()[0], 0.0);
  EXPECT_DOUBLE_EQ(moments[1].volumeMoments().centroid()[1], 0.0);
  EXPECT_DOUBLE_EQ(moments[1].volumeMoments().centroid()[2], 0.0);
  EXPECT_DOUBLE_EQ(moments[1].normal()[0], 1.0);
  EXPECT_DOUBLE_EQ(moments[1].normal()[1], 0.0);
  EXPECT_DOUBLE_EQ(moments[1].normal()[2], 0.0);
}

TEST(Partitioning, PartitionByNormalKMeansHarder) {
  ListedVolumeMoments<VolumeMomentsAndNormal> list;
  // Set distribution resulting from two shallow arcs. Symmetric across
  // y axis for simplicity.

  // Should be first partition
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(-0.2, 0.18, 0.0)),
                                 Normal(-0.2, 1.0, 0.0));
  // Should be second partition
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(0.2, -0.18, 0.0)),
                                 Normal(0.2, -1.0, 0.0));
  // Should be first partition
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(0.2, 0.18, 0.0)),
                                 Normal(0.2, 1.0, 0.0));
  // Should be second partition
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(-0.2, -0.18, 0.0)),
                                 Normal(-0.2, -1.0, 0.0));
  // Should be second partition
  list += VolumeMomentsAndNormal(VolumeMoments(0.1, Pt(0.0, -0.25, 0.0)),
                                 Normal(0.0, -1.0, 0.0));
  // Should be first partition
  list += VolumeMomentsAndNormal(VolumeMoments(0.1, Pt(0.0, 0.25, 0.0)),
                                 Normal(0.0, 1.0, 0.0));

  // Normalize all normals.
  for (auto& element : list) {
    element.normal().normalize();
  }
  list.multiplyByVolume();

  PartitionByNormal<ListedVolumeMoments<VolumeMomentsAndNormal>> partitioned(
      &list);
  partitioned.setup();
  auto iter = KMeans::partition(&partitioned);
  auto moments = partitioned.getPartitionedObjects();
  moments[0].normalizeByVolume();
  moments[1].normalizeByVolume();
  moments[0].normal().normalize();
  moments[1].normal().normalize();
  EXPECT_DOUBLE_EQ(moments[0].volumeMoments().volume(), 1.1);
  EXPECT_DOUBLE_EQ(moments[0].volumeMoments().centroid()[0], 0.0);
  EXPECT_NEAR(moments[0].volumeMoments().centroid()[1],
              (10 * 0.18 + 0.25) / 11.0, 1.0e-15);
  EXPECT_DOUBLE_EQ(moments[0].volumeMoments().centroid()[2], 0.0);
  EXPECT_DOUBLE_EQ(moments[0].normal()[0], 0.0);
  EXPECT_DOUBLE_EQ(moments[0].normal()[1], 1.0);
  EXPECT_DOUBLE_EQ(moments[0].normal()[2], 0.0);
  EXPECT_DOUBLE_EQ(moments[1].volumeMoments().volume(), 1.1);
  EXPECT_DOUBLE_EQ(moments[1].volumeMoments().centroid()[0], 0.0);
  EXPECT_NEAR(moments[1].volumeMoments().centroid()[1],
              -(10 * 0.18 + 0.25) / 11.0, 1.0e-15);
  EXPECT_DOUBLE_EQ(moments[1].volumeMoments().centroid()[2], 0.0);
  EXPECT_DOUBLE_EQ(moments[1].normal()[0], 0.0);
  EXPECT_DOUBLE_EQ(moments[1].normal()[1], -1.0);
  EXPECT_DOUBLE_EQ(moments[1].normal()[2], 0.0);
}

}  // namespace
