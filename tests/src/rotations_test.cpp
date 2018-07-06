// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/general/rotations.h"

#include <random>

#include "gtest/gtest.h"

#include "src/geometry/general/normal.h"
#include "src/helpers/helper.h"

namespace {

using namespace IRL;

TEST(Rotation, rotateNormalOntoNormal) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal normal_0 = Normal::normalized(random_normal(eng), random_normal(eng),
                                         random_normal(eng));
    Normal normal_1 = Normal::normalized(random_normal(eng), random_normal(eng),
                                         random_normal(eng));
    const UnitQuaternion rotation = rotateNormalOntoNormal(normal_0, normal_1);
    double rotation_amount;
    Normal rotation_axis;
    UnitQuaternion rotation_with_pointer = rotateNormalOntoNormal(
        normal_0, normal_1, &rotation_amount, &rotation_axis);
    const UnitQuaternion from_pointers(rotation_amount, rotation_axis);
    Normal rotated_normal = rotation * normal_0;
    EXPECT_EQ(rotated_normal == normal_1, true)
        << "Rotated normal: " << rotated_normal[0] << " " << rotated_normal[1]
        << " " << rotated_normal[2] << '\n'
        << "Target normal: " << normal_1[0] << " " << normal_1[1] << " "
        << normal_1[2] << '\n';
    EXPECT_DOUBLE_EQ(rotation[0], from_pointers[0]);
    EXPECT_DOUBLE_EQ(rotation[1], from_pointers[1]);
    EXPECT_DOUBLE_EQ(rotation[2], from_pointers[2]);
    EXPECT_DOUBLE_EQ(rotation[3], from_pointers[3]);
  }
}

TEST(Rotation, getOrthonormalSystem) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal normal = Normal::normalized(random_normal(eng), random_normal(eng),
                                       random_normal(eng));
    ReferenceFrame rotated_frame = getOrthonormalSystem(normal);
    Normal axis_0 = rotated_frame[0];
    Normal axis_1 = rotated_frame[1];
    Normal axis_2 = rotated_frame[2];
    EXPECT_NEAR(normal[0], axis_2[0], 1.0e-13);
    EXPECT_NEAR(normal[1], axis_2[1], 1.0e-13);
    EXPECT_NEAR(normal[2], axis_2[2], 1.0e-13);
    EXPECT_NEAR(axis_0 * axis_1, 0.0, 5.0e-16);
    EXPECT_NEAR(axis_0 * axis_2, 0.0, 5.0e-16);
    EXPECT_NEAR(axis_1 * axis_2, 0.0, 5.0e-16);
  }
}

TEST(Rotation, getSharedNormal) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal normal_0 = Normal::normalized(random_normal(eng), random_normal(eng),
                                         random_normal(eng));
    Normal normal_1 = Normal::normalized(random_normal(eng), random_normal(eng),
                                         random_normal(eng));
    double half_angle;
    Normal rotation_axis;
    Normal shared_normal =
        getSharedNormal(normal_1, normal_0, &half_angle, &rotation_axis);
    UnitQuaternion rotation_plus(half_angle, rotation_axis);
    UnitQuaternion rotation_minus(-half_angle, rotation_axis);
    Normal normal_plus = rotation_plus * shared_normal;
    Normal normal_minus = rotation_minus * shared_normal;
    EXPECT_NEAR(normal_plus[0], normal_0[0], 1.0e-13)
        << "Normal plus" << normal_plus[0] << " " << normal_plus[1] << " "
        << normal_plus[2] << '\n'
        << "Normal minus" << normal_minus[0] << " " << normal_minus[1] << " "
        << normal_minus[2] << '\n'
        << "Normal 0" << normal_0[0] << " " << normal_0[1] << " " << normal_0[2]
        << '\n'
        << "Normal 1" << normal_1[0] << " " << normal_1[1] << " " << normal_1[2]
        << "\n\n";
    EXPECT_NEAR(normal_plus[1], normal_0[1], 1.0e-13);
    EXPECT_NEAR(normal_plus[2], normal_0[2], 1.0e-13);
    EXPECT_NEAR(normal_minus[0], normal_1[0], 1.0e-13);
    EXPECT_NEAR(normal_minus[1], normal_1[1], 1.0e-13);
    EXPECT_NEAR(normal_minus[2], normal_1[2], 1.0e-13);
    Normal straight_shared_normal = getSharedNormal(normal_1, normal_0);
    EXPECT_EQ(straight_shared_normal == shared_normal, true);
  }

  Normal normal_0 = Normal::normalized(1.0, 0.0, 0.0);
  Normal normal_1 = Normal::normalized(1.0, 0.0, 0.0);
  double half_angle;
  Normal rotation_axis;
  Normal shared_normal =
      getSharedNormal(normal_1, normal_0, &half_angle, &rotation_axis);
  UnitQuaternion rotation_plus(half_angle, rotation_axis);
  UnitQuaternion rotation_minus(-half_angle, rotation_axis);
  Normal normal_plus = rotation_plus * shared_normal;
  Normal normal_minus = rotation_minus * shared_normal;
  EXPECT_NEAR(normal_plus[0], normal_0[0], 1.0e-13);
  EXPECT_NEAR(normal_plus[1], normal_0[1], 1.0e-13);
  EXPECT_NEAR(normal_plus[2], normal_0[2], 1.0e-13);
  EXPECT_NEAR(normal_minus[0], normal_1[0], 1.0e-13);
  EXPECT_NEAR(normal_minus[1], normal_1[1], 1.0e-13);
  EXPECT_NEAR(normal_minus[2], normal_1[2], 1.0e-13);

  normal_0 = Normal(0.0, 0.0, 1.0);
  normal_1 = Normal(0.0, 0.0, -1.0);
  shared_normal =
      getSharedNormal(normal_1, normal_0, &half_angle, &rotation_axis);
  rotation_plus = UnitQuaternion(half_angle, rotation_axis);
  rotation_minus = UnitQuaternion(-half_angle, rotation_axis);
  normal_plus = rotation_plus * shared_normal;
  normal_minus = rotation_minus * shared_normal;
  EXPECT_NEAR(normal_plus[0], normal_0[0], 1.0e-13);
  EXPECT_NEAR(normal_plus[1], normal_0[1], 1.0e-13);
  EXPECT_NEAR(normal_plus[2], normal_0[2], 1.0e-13);
  EXPECT_NEAR(normal_minus[0], normal_1[0], 1.0e-13);
  EXPECT_NEAR(normal_minus[1], normal_1[1], 1.0e-13);
  EXPECT_NEAR(normal_minus[2], normal_1[2], 1.0e-13);
}

}  // namespace
