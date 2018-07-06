// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/helpers/geometric_cutting_helpers.h"

#include <float.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include "src/planar_reconstruction/planar_localizer.h"
#include "src/planar_reconstruction/planar_reconstruction.h"
#include "src/planar_reconstruction/planar_separator.h"

#include "gtest/gtest.h"
#include "src/generic_cutting/recursive_simplex_cutting/lookup_tables.h"

#include "src/geometry/general/plane.h"
#include "src/geometry/polygons/polygon.h"
#include "src/helpers/helper.h"
#include "src/helpers/mymath.h"

namespace {

using namespace IRL;

TEST(GeometricCuttingHelpers, distanceBetweenPts) {
  // Check correct calculation of distance between two given points
  Pt pt0(2.0, 2.0, 2.0);
  Pt pt1(1.0, 1.0, 1.0);
  EXPECT_DOUBLE_EQ(distanceBetweenPts(pt0, pt1), std::sqrt(3.0));
  EXPECT_DOUBLE_EQ(distanceBetweenPts(pt1, pt0), std::sqrt(3.0));
}

TEST(GeometricCuttingHelpers, squaredDistanceBetweenPts) {
  // Check correct calculation of distance between two given points
  Pt pt0(2.0, 2.0, 2.0);
  Pt pt1(1.0, 1.0, 1.0);
  EXPECT_DOUBLE_EQ(squaredDistanceBetweenPts(pt0, pt1), 3.0);
  EXPECT_DOUBLE_EQ(squaredDistanceBetweenPts(pt1, pt0), 3.0);
}

TEST(GeometricCutting, isPtInternal) {
  PlanarLocalizer localizer_sheet;
  localizer_sheet.addPlane(Plane(Normal(0.0, -1.0, 0.0), 0.25));
  localizer_sheet.addPlane(Plane(Normal(0.0, 1.0, 0.0), 0.25));
  EXPECT_EQ(isPtInternal(Pt(0.25, -0.5, 0.8), localizer_sheet), false);
  EXPECT_EQ(isPtInternal(Pt(12.0, 0.0, -0.8), localizer_sheet), true);
  EXPECT_EQ(isPtInternal(Pt(62.0, 0.5, 0.3), localizer_sheet), false);

  PlanarSeparator inverse_sheet;
  inverse_sheet.addPlane(Plane(Normal(0.0, -1.0, 0.0), -0.25));
  inverse_sheet.addPlane(Plane(Normal(0.0, 1.0, 0.0), -0.25));
  inverse_sheet.flipCutting();
  EXPECT_EQ(isPtInternal(Pt(0.25, -0.5, 0.8), inverse_sheet), true);
  EXPECT_EQ(isPtInternal(Pt(12.0, 0.0, -0.8), inverse_sheet), false);
  EXPECT_EQ(isPtInternal(Pt(62.0, 0.5, 0.3), inverse_sheet), true);

  PlanarSeparator inverse_triangle;
  inverse_triangle.addPlane(
      Plane(Normal(-0.5 * std::sqrt(2.0), 0.5 * std::sqrt(2.0), 0.0), 0.0));
  inverse_triangle.addPlane(
      Plane(Normal(0.5 * std::sqrt(2.0), 0.5 * std::sqrt(2.0), 0.0), 0.0));
  inverse_triangle.flipCutting();
  EXPECT_EQ(isPtInternal(Pt(0.25, -0.5, 0.8), inverse_triangle), true);
  EXPECT_EQ(isPtInternal(Pt(0.0, -1.0e-14, 0.0), inverse_triangle), true);
  EXPECT_EQ(isPtInternal(Pt(-50.0, 0.5, -0.8), inverse_triangle), true);
  EXPECT_EQ(isPtInternal(Pt(50.0, 0.5, -0.8), inverse_triangle), true);
  EXPECT_EQ(isPtInternal(Pt(0.0, 0.5, 0.3), inverse_triangle), false);
}

TEST(GeometricCuttingHelpers, getGeometricCaseId) {
  // Just test a handful of cases to make sure correct case is identified from
  // signs
  std::array<double, 4> distance = {-1.0, -1.0, -1.0, -1.0};
  int my_case = 0;
  EXPECT_EQ(getGeometricCaseId(distance), 0);
  for (UnsignedIndex_t v = 0; v < 4; ++v) {
    distance[v] = -distance[v];
    my_case = my_case + static_cast<int>(std::pow(2, v));
    EXPECT_EQ(getGeometricCaseId(distance), my_case);
  }
  EXPECT_EQ(getGeometricCaseId(distance), 15);
  distance[2] = -1.0;
  EXPECT_EQ(getGeometricCaseId(distance), 11);
  distance[3] = -1.0;
  EXPECT_EQ(getGeometricCaseId(distance), 3);
}

}  // namespace
