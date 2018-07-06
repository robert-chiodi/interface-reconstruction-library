// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/polyhedrons/rectangular_cuboid.h"

#include "gtest/gtest.h"

#include <algorithm>

#include "src/geometry/general/pt.h"
#include "src/geometry/general/rotations.h"
#include "src/geometry/polygons/tri.h"
#include "src/helpers/geometric_cutting_helpers.h"

namespace {

using namespace IRL;

TEST(Polyhedra, RectangularCuboid) {
  // Construction and access methods
  RectangularCuboid two_point_set_rect_cubic =
      RectangularCuboid::fromBoundingPts(Pt(-0.5, -0.5, -0.5),
                                         Pt(0.5, 0.5, 0.5));

  double flat_pt_list[24] = {0.5, -0.5, -0.5, 0.5, 0.5,  -0.5, 0.5,  0.5,
                             0.5, 0.5,  -0.5, 0.5, -0.5, -0.5, -0.5, -0.5,
                             0.5, -0.5, -0.5, 0.5, 0.5,  -0.5, -0.5, 0.5};

  RectangularCuboid flat_list_cuboid =
      RectangularCuboid::fromRawDoublePointer(8, flat_pt_list);
  for (UnsignedIndex_t v = 0; v < 8; ++v) {
    EXPECT_DOUBLE_EQ(flat_list_cuboid[v].x(), two_point_set_rect_cubic[v].x());
    EXPECT_DOUBLE_EQ(flat_list_cuboid[v].y(), two_point_set_rect_cubic[v].y());
    EXPECT_DOUBLE_EQ(flat_list_cuboid[v].z(), two_point_set_rect_cubic[v].z());
  }

  // Shifting cubes
  auto rect_cubic = two_point_set_rect_cubic;
  rect_cubic.shift(-1.0, -1.0, -1.0);
  for (UnsignedIndex_t v = 0; v < 8; ++v) {
    EXPECT_DOUBLE_EQ(rect_cubic[v].x(), two_point_set_rect_cubic[v].x() - 1.0);
    EXPECT_DOUBLE_EQ(rect_cubic[v].y(), two_point_set_rect_cubic[v].y() - 1.0);
    EXPECT_DOUBLE_EQ(rect_cubic[v].z(), two_point_set_rect_cubic[v].z() - 1.0);
  }
  EXPECT_DOUBLE_EQ(rect_cubic.calculateSideLength(0), 1.0);
  EXPECT_DOUBLE_EQ(rect_cubic.calculateSideLength(1), 1.0);
  EXPECT_DOUBLE_EQ(rect_cubic.calculateSideLength(2), 1.0);

  rect_cubic.shift(1.0, 1.0, 1.0);
  for (UnsignedIndex_t v = 0; v < 8; ++v) {
    EXPECT_DOUBLE_EQ(rect_cubic[v].x(), two_point_set_rect_cubic[v].x());
    EXPECT_DOUBLE_EQ(rect_cubic[v].y(), two_point_set_rect_cubic[v].y());
    EXPECT_DOUBLE_EQ(rect_cubic[v].z(), two_point_set_rect_cubic[v].z());
  }

  // Calculation of geomtric quantities
  EXPECT_DOUBLE_EQ(two_point_set_rect_cubic.calculateVolume(), 1.0);
  Pt centroid = two_point_set_rect_cubic.calculateCentroid();
  EXPECT_DOUBLE_EQ(centroid.x(), 0.0);
  EXPECT_DOUBLE_EQ(centroid.y(), 0.0);
  EXPECT_DOUBLE_EQ(centroid.z(), 0.0);
}
}  // namespace
