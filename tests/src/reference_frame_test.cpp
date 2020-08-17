// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/general/reference_frame.h"

#include <random>

#include "gtest/gtest.h"

#include "src/geometry/general/normal.h"
#include "src/helpers/helper.h"

namespace {

using namespace IRL;

TEST(Rotation, ReferenceFrame) {
  // Check counstruction and access
  Normal n1 = Normal::normalized(1.0, 0.0, 0.0);
  Normal n2 = Normal::normalized(0.0, 1.0, 0.0);
  Normal n3 = Normal::normalized(0.0, 0.0, 1.0);
  ReferenceFrame a_ref_frame(n1, n2, n3);
  Normal take_from_frame = a_ref_frame[0];
  EXPECT_NEAR(take_from_frame[0], 1.0, DBL_EPSILON);
  EXPECT_NEAR(take_from_frame[1], 0.0, DBL_EPSILON);
  EXPECT_NEAR(take_from_frame[2], 0.0, DBL_EPSILON);
  take_from_frame = a_ref_frame[1];
  EXPECT_NEAR(take_from_frame[0], 0.0, DBL_EPSILON);
  EXPECT_NEAR(take_from_frame[1], 1.0, DBL_EPSILON);
  EXPECT_NEAR(take_from_frame[2], 0.0, DBL_EPSILON);
  a_ref_frame[1] = a_ref_frame[2];
  take_from_frame = a_ref_frame[1];
  EXPECT_NEAR(take_from_frame[0], 0.0, DBL_EPSILON);
  EXPECT_NEAR(take_from_frame[1], 0.0, DBL_EPSILON);
  EXPECT_NEAR(take_from_frame[2], 1.0, DBL_EPSILON);
}
}  // namespace
