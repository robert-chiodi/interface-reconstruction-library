// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/general/normal.h"

#include <algorithm>
#include <random>

#include "gtest/gtest.h"

#include "src/geometry/general/rotations.h"

#include "src/helpers/helper.h"
namespace {

using namespace IRL;

TEST(GeometryClasses, Normal) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.

  // Check constructing, setting, and normalizing
  static const int ncycles = 10;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);

  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal n1 = Normal::normalized(random_normal(eng), random_normal(eng),
                                   random_normal(eng));
    EXPECT_DOUBLE_EQ(n1 * n1, 1.0);
  }
  Pt pt1(-1.1, 1.1, -1.1);
  Pt pt2(-1.0, 1.0, -1.0);
  Normal ptdiff = Normal::fromPtNormalized(pt1 - pt2);
  EXPECT_DOUBLE_EQ(ptdiff[0], -std::sqrt(3.0) / 3.0);
  EXPECT_DOUBLE_EQ(ptdiff[1], std::sqrt(3.0) / 3.0);
  EXPECT_DOUBLE_EQ(ptdiff[2], -std::sqrt(3.0) / 3.0);

  Normal n1;
  n1 = Normal(0.0, 1.0, 1.0);
  n1.normalize();
  EXPECT_DOUBLE_EQ(n1[0], 0.0);
  EXPECT_DOUBLE_EQ(n1[1], 0.5 * std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(n1[2], 0.5 * std::sqrt(2.0));
  n1 = -n1;
  EXPECT_DOUBLE_EQ(n1[0], 0.0);
  EXPECT_DOUBLE_EQ(n1[1], -0.5 * std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(n1[2], -0.5 * std::sqrt(2.0));

  n1 = Normal(1.0, 0.0, -1.0);
  n1.normalize();
  EXPECT_DOUBLE_EQ(n1[0], 0.5 * std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(n1[1], 0.0);
  EXPECT_DOUBLE_EQ(n1[2], -0.5 * std::sqrt(2.0));

  // Check dot products with points
  Pt pt0(0.25, 0.3, -0.25);
  EXPECT_DOUBLE_EQ(pt0 * n1, 2 * 0.25 * 0.5 * std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(n1 * pt0, 2 * 0.25 * 0.5 * std::sqrt(2.0));

  // Check dot product with normals
  double input_normal[3] = {-1.0, 0.0, 1.0};
  Normal n2 = Normal::fromRawDoublePointerNormalized(input_normal);
  EXPECT_DOUBLE_EQ(n1 * n2, -1.0);
  EXPECT_DOUBLE_EQ(n2 * n1, -1.0);

  // Check overload of /=
  n2 /= 0.5 * std::sqrt(2.0);
  EXPECT_DOUBLE_EQ(n2[0], -1.0);
  EXPECT_DOUBLE_EQ(n2[1], 0.0);
  EXPECT_DOUBLE_EQ(n2[2], 1.0);

  // Check overload of double*Normal = Normal
  n1 = Normal(1.0, -1.0, 1.0);
  n1.normalize();
  Normal mult_normal = 0.5 * n1;
  EXPECT_DOUBLE_EQ(mult_normal[0], std::sqrt(3.0) / 6.0);
  EXPECT_DOUBLE_EQ(mult_normal[1], -std::sqrt(3.0) / 6.0);
  EXPECT_DOUBLE_EQ(mult_normal[2], std::sqrt(3.0) / 6.0);
  mult_normal = n1 * 0.5;
  EXPECT_DOUBLE_EQ(mult_normal[0], std::sqrt(3.0) / 6.0);
  EXPECT_DOUBLE_EQ(mult_normal[1], -std::sqrt(3.0) / 6.0);
  EXPECT_DOUBLE_EQ(mult_normal[2], std::sqrt(3.0) / 6.0);

  // Test boolean
  n1 = Normal(1.0, 0.0, 0.0);
  n2 = Normal(1.0, 0.0, 0.0);
  EXPECT_EQ(n1 == n2, true);
  n2 = Normal(0.0, 1.0, 0.0);
  EXPECT_EQ(n1 == n2, false);
  double angle = std::acos(global_constants::SAME_VEC) +
                 1.0e-4 * global_constants::SAME_VEC;
  UnitQuaternion rotation(angle, Normal(0.0, 0.0, 1.0));
  n2 = rotation * n1;
  EXPECT_EQ(n1 == n2, false);
  EXPECT_EQ(n1 != n2, true);
  angle = std::acos(global_constants::SAME_VEC) -
          1.0e-4 * global_constants::SAME_VEC;
  UnitQuaternion rotation2(angle, Normal(0.0, 0.0, 1.0));
  n2 = rotation2 * n1;
  EXPECT_EQ(n1 == n2, true);
  EXPECT_EQ(n1 != n2, false);

  // Test serialization
  ByteBuffer buffer;
  std::array<Normal, 2> normal_to_pack;
  normal_to_pack[0] = Normal::normalized(2.8, -4.3, 8.2);
  normal_to_pack[1] = Normal::normalized(6.9, 2.2, -4.03);
  serializeAndPack(normal_to_pack[0], &buffer);
  serializeAndPack(normal_to_pack[1], &buffer);
  buffer.resetBufferPointer();
  std::array<Normal, 2> recv_normal;
  unpackAndStore(&recv_normal[0], &buffer);
  unpackAndStore(&recv_normal[1], &buffer);
  for (UnsignedIndex_t n = 0; n < normal_to_pack.size(); ++n) {
    EXPECT_DOUBLE_EQ(recv_normal[n][0], normal_to_pack[n][0]);
    EXPECT_DOUBLE_EQ(recv_normal[n][1], normal_to_pack[n][1]);
    EXPECT_DOUBLE_EQ(recv_normal[n][2], normal_to_pack[n][2]);
  }
}
}  // namespace
