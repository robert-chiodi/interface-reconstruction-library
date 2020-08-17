// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/general/plane.h"

#include <algorithm>
#include <random>

#include "gtest/gtest.h"

#include "src/geometry/general/rotations.h"

#include "src/helpers/helper.h"
namespace {

using namespace IRL;

TEST(GeometryClasses, Plane) {
  // Setting and accessing
  Normal n1 = Normal::normalized(1.0, 1.0, 0.0);
  Plane plane_1(n1, 0.5);
  EXPECT_DOUBLE_EQ(plane_1.normal()[0], 0.5 * std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(plane_1.normal()[1], 0.5 * std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(plane_1.normal()[2], 0.0);
  EXPECT_DOUBLE_EQ(plane_1.distance(), 0.5);
  plane_1.distance() = 1.0;
  EXPECT_DOUBLE_EQ(plane_1.distance(), 1.0);
  Normal out_norm = plane_1.normal();
  EXPECT_DOUBLE_EQ(out_norm[0], 0.5 * std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(out_norm[1], 0.5 * std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(out_norm[2], 0.0);

  double init_normal[3] = {1.0, 1.0, 0.0};
  Plane plane_2(Normal::fromRawDoublePointerNormalized(init_normal), 0.5);
  EXPECT_DOUBLE_EQ(plane_1.normal()[0], plane_2.normal()[0]);
  EXPECT_DOUBLE_EQ(plane_1.normal()[1], plane_2.normal()[1]);
  EXPECT_DOUBLE_EQ(plane_1.normal()[2], plane_2.normal()[2]);

  // Reset normal for a plane
  Normal new_normal = Normal::normalized(0.0, 1.0, 1.0);
  plane_1.normal() = new_normal;
  plane_1.distance() = 8.4;
  EXPECT_DOUBLE_EQ(plane_1.normal()[0], 0.0);
  EXPECT_DOUBLE_EQ(plane_1.normal()[1], 0.5 * std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(plane_1.normal()[2], 0.5 * std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(plane_1.distance(), 8.4);
  plane_1.normal() = Normal::normalized(1.0, 0.0, 1.0);
  EXPECT_DOUBLE_EQ(plane_1.normal()[0], 0.5 * std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(plane_1.normal()[1], 0.0);
  EXPECT_DOUBLE_EQ(plane_1.normal()[2], 0.5 * std::sqrt(2.0));
  plane_1.normal() = Normal::fromRawDoublePointerNormalized(init_normal);
  EXPECT_DOUBLE_EQ(plane_1.normal()[0], 0.5 * std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(plane_1.normal()[1], 0.5 * std::sqrt(2.0));

  // Getting signed distance that a point is from a plane
  Pt pt0(0.0, 1.0, 0.0);
  Plane cut_plane(Normal(0.0, 1.0, 0.0), 0.5);
  EXPECT_DOUBLE_EQ(cut_plane.signedDistanceToPoint(pt0), 0.5);
  pt0[1] = -pt0[1];
  EXPECT_DOUBLE_EQ(cut_plane.signedDistanceToPoint(pt0), -1.5);
  Pt pt1(1.0, 0.0, 0.0);
  EXPECT_DOUBLE_EQ(cut_plane.signedDistanceToPoint(pt1), -0.5);
  Plane cut_plane_2(Normal(-1.0, 0.0, 0.0), -0.5);
  EXPECT_DOUBLE_EQ(cut_plane_2.signedDistanceToPoint(pt1), -0.5);
  Plane cut_plane_3(Normal(1.0, 0.0, 0.0), -0.5);
  EXPECT_DOUBLE_EQ(cut_plane_3.signedDistanceToPoint(pt1), 1.5);

  // Test serialization
  ByteBuffer buffer;
  std::array<Plane, 2> plane_to_pack;
  plane_to_pack[0] = Plane(Normal::normalized(2.8, -4.3, 8.2), -42.8);
  plane_to_pack[1] = Plane(Normal::normalized(6.9, 2.2, -4.03), 20.8);
  serializeAndPack(plane_to_pack[0], &buffer);
  serializeAndPack(plane_to_pack[1], &buffer);
  buffer.resetBufferPointer();
  std::array<Plane, 2> recv_plane;
  unpackAndStore(&recv_plane[0], &buffer);
  unpackAndStore(&recv_plane[1], &buffer);
  for (UnsignedIndex_t n = 0; n < plane_to_pack.size(); ++n) {
    EXPECT_DOUBLE_EQ(recv_plane[n].normal()[0], plane_to_pack[n].normal()[0]);
    EXPECT_DOUBLE_EQ(recv_plane[n].normal()[1], plane_to_pack[n].normal()[1]);
    EXPECT_DOUBLE_EQ(recv_plane[n].normal()[2], plane_to_pack[n].normal()[2]);
    EXPECT_DOUBLE_EQ(recv_plane[n].distance(), plane_to_pack[n].distance());
  }
}

}  // namespace
