// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/generic_cutting/generic_cutting.h"

#include "gtest/gtest.h"

#include <algorithm>

#include "src/geometry/polyhedrons/rectangular_cuboid.h"

namespace {

using namespace IRL;

TEST(SimplexCutting, SplitTest) {
  RectangularCuboid volume_to_cut = RectangularCuboid::fromBoundingPts(
      Pt(-10.5, -10.5, -10.5), Pt(-10.0, -10.0, -10.0));
  PlanarSeparator separator =
      PlanarSeparator::fromOnePlane(Plane(Normal(1.0, 0.0, 0.0), 0.0));
  auto volume_moments =
      getNormalizedVolumeMoments<VolumeMoments, SimplexCutting>(volume_to_cut,
                                                                separator);
  EXPECT_NEAR(volume_moments.volume(), 0.125, 1.0e-15);
  EXPECT_NEAR(volume_moments.centroid()[0], -10.25, 1.0e-15);
  EXPECT_NEAR(volume_moments.centroid()[1], -10.25, 1.0e-15);
  EXPECT_NEAR(volume_moments.centroid()[2], -10.25, 1.0e-15);
}

}  // namespace
