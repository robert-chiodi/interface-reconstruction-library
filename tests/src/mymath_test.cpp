// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/helpers/mymath.h"

#include <float.h>
#include <cmath>

#include "gtest/gtest.h"

#include "src/geometry/general/pt.h"

namespace {

using namespace IRL;

TEST(MyMath, deg2Rad) {
  EXPECT_DOUBLE_EQ(deg2Rad(45.0), 0.25 * M_PI);
  EXPECT_DOUBLE_EQ(deg2Rad(-45.0), -0.25 * M_PI);
  EXPECT_DOUBLE_EQ(deg2Rad(0.0), 0.0);
  EXPECT_DOUBLE_EQ(deg2Rad(360.0), 2.0 * M_PI);
  EXPECT_DOUBLE_EQ(deg2Rad(135.0), 0.75 * M_PI);
}

TEST(MyMath, rad2Deg) {
  EXPECT_DOUBLE_EQ(rad2Deg(0.25 * M_PI), 45.0);
  EXPECT_DOUBLE_EQ(rad2Deg(-0.25 * M_PI), -45.0);
  EXPECT_DOUBLE_EQ(rad2Deg(0.0), 0.0);
  EXPECT_DOUBLE_EQ(rad2Deg(2.0 * M_PI), 360.0);
  EXPECT_DOUBLE_EQ(rad2Deg(0.75 * M_PI), 135.0);
}

TEST(MyMath, angleNormalize) {
  EXPECT_NEAR(angleNormalize(deg2Rad(540)), deg2Rad(180), 5.0 * DBL_EPSILON);
  EXPECT_NEAR(angleNormalize(deg2Rad(-10)), deg2Rad(350), 5.0 * DBL_EPSILON);
  EXPECT_NEAR(angleNormalize(deg2Rad(-90)), deg2Rad(270), 5.0 * DBL_EPSILON);
  EXPECT_NEAR(angleNormalize(deg2Rad(720)), deg2Rad(0.0), 5.0 * DBL_EPSILON);
}

TEST(MyMath, signedAngleNormalize) {
  EXPECT_NEAR(signedAngleNormalize(deg2Rad(540)), deg2Rad(180),
              5.0 * DBL_EPSILON);
  EXPECT_NEAR(signedAngleNormalize(deg2Rad(-10)), deg2Rad(-10),
              5.0 * DBL_EPSILON);
  EXPECT_NEAR(signedAngleNormalize(deg2Rad(-90)), deg2Rad(-90),
              5.0 * DBL_EPSILON);
  EXPECT_NEAR(signedAngleNormalize(deg2Rad(720)), deg2Rad(0.0),
              5.0 * DBL_EPSILON);
  EXPECT_NEAR(signedAngleNormalize(deg2Rad(-540)), deg2Rad(-180.0),
              5.0 * DBL_EPSILON);
}

TEST(MyMath, crossProduct) {
  Pt pt1(1.0, 0.0, 0.0);
  Pt pt2(0.0, 1.0, 0.0);
  Pt result = crossProduct(pt1, pt2);
  EXPECT_DOUBLE_EQ(magnitude(result), 1.0);
  EXPECT_DOUBLE_EQ(result[0], 0.0);
  EXPECT_DOUBLE_EQ(result[1], 0.0);
  EXPECT_DOUBLE_EQ(result[2], 1.0);
  result = crossProduct(pt2, pt1);
  EXPECT_DOUBLE_EQ(magnitude(result), 1.0);
  EXPECT_DOUBLE_EQ(result[0], 0.0);
  EXPECT_DOUBLE_EQ(result[1], 0.0);
  EXPECT_DOUBLE_EQ(result[2], -1.0);
  pt1 = Pt(0.0, 0.0, 1.0);
  pt2 = Pt(1.0, 0.0, 0.0);
  result = crossProduct(pt1, pt2);
  EXPECT_DOUBLE_EQ(magnitude(result), 1.0);
  EXPECT_DOUBLE_EQ(result[0], 0.0);
  EXPECT_DOUBLE_EQ(result[1], 1.0);
  EXPECT_DOUBLE_EQ(result[2], 0.0);
}

TEST(MyMath, crossProductNormalized) {
  Pt pt1(0.5, 0.0, 0.0);
  Pt pt2(0.0, 0.5, 0.0);
  Pt result = crossProductNormalized(pt1, pt2);
  EXPECT_DOUBLE_EQ(magnitude(result), 1.0);
  EXPECT_DOUBLE_EQ(result[0], 0.0);
  EXPECT_DOUBLE_EQ(result[1], 0.0);
  EXPECT_DOUBLE_EQ(result[2], 1.0);
  result = crossProductNormalized(pt2, pt1);
  EXPECT_DOUBLE_EQ(magnitude(result), 1.0);
  EXPECT_DOUBLE_EQ(result[0], 0.0);
  EXPECT_DOUBLE_EQ(result[1], 0.0);
  EXPECT_DOUBLE_EQ(result[2], -1.0);
  pt1 = Pt(0.0, 0.0, 0.5);
  pt2 = Pt(0.5, 0.0, 0.0);
  result = crossProductNormalized(pt1, pt2);
  EXPECT_DOUBLE_EQ(magnitude(result), 1.0);
  EXPECT_DOUBLE_EQ(result[0], 0.0);
  EXPECT_DOUBLE_EQ(result[1], 1.0);
  EXPECT_DOUBLE_EQ(result[2], 0.0);
}

TEST(MyMath, magnitude) {
  Pt pt0(0.0, 0.0, 0.0);
  Pt pt1(1.0, 0.0, 0.0);
  Pt pt2(1.0, 1.0, 0.0);
  EXPECT_DOUBLE_EQ(magnitude(pt2 - pt0), std::sqrt(2.0));
  EXPECT_DOUBLE_EQ(magnitude(pt2 - pt1), 1.0);
  EXPECT_DOUBLE_EQ(magnitude(pt1 - pt0), 1.0);
}

TEST(MyMath, squaredMagnitude) {
  Pt pt0(0.0, 0.0, 0.0);
  Pt pt1(1.0, 0.0, 0.0);
  Pt pt2(1.0, 1.0, 0.0);
  EXPECT_DOUBLE_EQ(squaredMagnitude(pt2 - pt0), 2.0);
  EXPECT_DOUBLE_EQ(squaredMagnitude(pt2 - pt1), 1.0);
  EXPECT_DOUBLE_EQ(squaredMagnitude(pt1 - pt0), 1.0);
}

}  // namespace
