// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/general/unit_quaternion.h"

#include <random>

#include "gtest/gtest.h"

#include "src/geometry/general/normal.h"
#include "src/helpers/helper.h"

namespace {

using namespace IRL;

TEST(Rotation, UnitQuaternion) {
  Normal x_normal = Normal::normalized(1.0, 0.0, 0.0);
  Normal z_rotation = Normal::normalized(0.0, 0.0, 1.0);
  // Test rotation of normals
  double deg = {10.0};
  double deg_sum = {0.0};
  UnitQuaternion rotate_z(deg2Rad(deg), z_rotation);
  UnitQuaternion summed_rotation_xy(0.0, z_rotation);
  Normal rotated_normal;
  for (int i = 0; i < 36; ++i) {
    deg_sum += deg;
    summed_rotation_xy = rotate_z * summed_rotation_xy;
    rotated_normal = summed_rotation_xy * x_normal;
    EXPECT_NEAR(rotated_normal[0], std::cos(deg2Rad(deg_sum)),
                10.0 * DBL_EPSILON);
    EXPECT_NEAR(rotated_normal[1], std::sin(deg2Rad(deg_sum)),
                10.0 * DBL_EPSILON);
    EXPECT_NEAR(rotated_normal[2], 0.0, DBL_EPSILON);
  }
}
}  // namespace
