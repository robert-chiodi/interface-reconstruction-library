// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/general/pt.h"

#include <algorithm>
#include <random>

#include "gtest/gtest.h"

#include "src/geometry/general/rotations.h"

#include "src/geometry/general/plane.h"
#include "src/helpers/helper.h"
namespace {

using namespace IRL;

TEST(GeometryClasses, Pt) {
  // Construction and access
  Pt pt0;
  pt0 = Pt(2.0, 1.0, 0.0);
  Pt pt1(0.0, 1.0, 2.0);
  double init[3] = {2.0, 1.0, 0.0};
  Pt array_pt = Pt::fromRawDoublePointer(init);
  EXPECT_DOUBLE_EQ(pt0[0], array_pt[0]);
  EXPECT_DOUBLE_EQ(pt0[1], array_pt[1]);
  EXPECT_DOUBLE_EQ(pt0[2], array_pt[2]);

  // Operators
  Pt diff_pt = pt0 - pt1;
  EXPECT_DOUBLE_EQ(diff_pt.x(), 2.0);
  EXPECT_DOUBLE_EQ(diff_pt.y(), 0.0);
  EXPECT_DOUBLE_EQ(diff_pt.z(), -2.0);

  Pt plus_pt = pt0 + pt1;
  EXPECT_DOUBLE_EQ(plus_pt.x(), 2.0);
  EXPECT_DOUBLE_EQ(plus_pt.y(), 2.0);
  EXPECT_DOUBLE_EQ(plus_pt.z(), 2.0);

  Pt right_mult_pt = pt0 * 4.0;
  EXPECT_DOUBLE_EQ(right_mult_pt.x(), 8.0);
  EXPECT_DOUBLE_EQ(right_mult_pt.y(), 4.0);
  EXPECT_DOUBLE_EQ(right_mult_pt.z(), 0.0);

  Pt left_mult_pt = 4.0 * pt0;
  EXPECT_DOUBLE_EQ(right_mult_pt.x(), left_mult_pt.x());
  EXPECT_DOUBLE_EQ(right_mult_pt.y(), left_mult_pt.y());
  EXPECT_DOUBLE_EQ(right_mult_pt.z(), left_mult_pt.z());

  Pt divide_pt = pt0 / 2.0;
  EXPECT_DOUBLE_EQ(divide_pt.x(), 1.0);
  EXPECT_DOUBLE_EQ(divide_pt.y(), 0.5);
  EXPECT_DOUBLE_EQ(divide_pt.z(), 0.0);

  // Test serialization
  ByteBuffer buffer;
  std::array<Pt, 2> pt_to_pack;
  pt_to_pack[0] = Pt(2.8, -4.3, 8.2);
  pt_to_pack[1] = Pt(6.9, 2.2, -50.03);
  serializeAndPack(pt_to_pack[0], &buffer);
  serializeAndPack(pt_to_pack[1], &buffer);
  buffer.resetBufferPointer();
  std::array<Pt, 2> recv_pt;
  unpackAndStore(&recv_pt[0], &buffer);
  unpackAndStore(&recv_pt[1], &buffer);
  for (UnsignedIndex_t n = 0; n < pt_to_pack.size(); ++n) {
    EXPECT_DOUBLE_EQ(recv_pt[n].x(), pt_to_pack[n].x());
    EXPECT_DOUBLE_EQ(recv_pt[n].y(), pt_to_pack[n].y());
    EXPECT_DOUBLE_EQ(recv_pt[n].z(), pt_to_pack[n].z());
  }
}

TEST(Pt, ptPlaneIntersectsLine) {
  // Check the correct location of the point where a plane intersects a line
  Pt pt0(0.0, 0.0, 0.0);
  Pt pt1(1.0, 0.0, 0.0);
  Plane easy_cut_plane(Normal(1.0, 0.0, 0.0), 0.5);
  Pt easy_intersection_pt =
      Pt::fromEdgeIntersection(pt0, easy_cut_plane.signedDistanceToPoint(pt0),
                               pt1, easy_cut_plane.signedDistanceToPoint(pt1));
  EXPECT_DOUBLE_EQ(easy_intersection_pt.x(), 0.5);
  EXPECT_DOUBLE_EQ(easy_intersection_pt.y(), 0.0);
  EXPECT_DOUBLE_EQ(easy_intersection_pt.z(), 0.0);

  Pt pt2(0.0, 0.0, 0.0);
  Pt pt3(1.0, 1.0, 0.5);
  Plane cut_plane(Normal(0.5 * std::sqrt(2.0), 0.5 * std::sqrt(2.0), 0.0), 0.5);
  Pt intersection_pt =
      Pt::fromEdgeIntersection(pt2, cut_plane.signedDistanceToPoint(pt2), pt3,
                               cut_plane.signedDistanceToPoint(pt3));
  EXPECT_DOUBLE_EQ(intersection_pt.x(),
                   0.5 / (0.5 * std::sqrt(2) + 0.5 * std::sqrt(2.0)));
  EXPECT_DOUBLE_EQ(intersection_pt.y(),
                   0.5 / (0.5 * std::sqrt(2) + 0.5 * std::sqrt(2.0)));
  EXPECT_DOUBLE_EQ(intersection_pt.z(),
                   0.25 / (0.5 * std::sqrt(2) + 0.5 * std::sqrt(2.0)));
}

}  // namespace
