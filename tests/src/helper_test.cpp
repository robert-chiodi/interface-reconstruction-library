// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/helpers/helper.h"

#include <float.h>
#include <algorithm>
#include <array>
#include <cmath>

#include "gtest/gtest.h"

namespace {

using namespace IRL;

TEST(Helper, safelyTiny) {
  double val = -DBL_MIN * 0.1;
  EXPECT_DOUBLE_EQ(safelyTiny(val), -DBL_MIN);
  val = -val;
  EXPECT_DOUBLE_EQ(safelyTiny(val), DBL_MIN);
  val = -14.0;
  EXPECT_DOUBLE_EQ(safelyTiny(val), val);
  val = -val;
  EXPECT_DOUBLE_EQ(safelyTiny(val), val);
}

TEST(Helper, safelyEpsilon) {
  double val = -DBL_EPSILON * 0.1;
  EXPECT_DOUBLE_EQ(safelyEpsilon(val), -DBL_EPSILON);
  val = -val;
  EXPECT_DOUBLE_EQ(safelyEpsilon(val), DBL_EPSILON);
  val = -14.0;
  EXPECT_DOUBLE_EQ(safelyEpsilon(val), val);
  val = -val;
  EXPECT_DOUBLE_EQ(safelyEpsilon(val), val);
}

TEST(Helper, safelySmall) {
  double val = -0.8;
  EXPECT_DOUBLE_EQ(safelySmall(val, 0.4), val);
  EXPECT_DOUBLE_EQ(safelySmall(val, 1.0), -1.0);
  val = -val;
  EXPECT_DOUBLE_EQ(safelySmall(val, 0.4), val);
  EXPECT_DOUBLE_EQ(safelySmall(val, 1.0), 1.0);
}

TEST(Helper, clipBetween) {
  double val = 25.0;
  EXPECT_DOUBLE_EQ(clipBetween(0.0, val, 100.0), val);
  EXPECT_DOUBLE_EQ(clipBetween(50.0, val, 100.0), 50.0);
  EXPECT_DOUBLE_EQ(clipBetween(0.0, val, 20.0), 20);
  EXPECT_DOUBLE_EQ(clipBetween(-80.0, val, 100.0), val);
  val = -50.0;
  EXPECT_DOUBLE_EQ(clipBetween(-100.0, val, 0.0), val);
  EXPECT_DOUBLE_EQ(clipBetween(-40.0, val, 0.0), -40.0);
  EXPECT_DOUBLE_EQ(clipBetween(-100.0, val, -60.0), -60.0);
}

TEST(Helper, wantPurelyInternal) {
  EXPECT_EQ(wantPurelyInternal(0.0), false);
  EXPECT_EQ(wantPurelyInternal(0.4), false);
  EXPECT_EQ(wantPurelyInternal(global_constants::VF_HIGH - DBL_EPSILON), false);
  EXPECT_EQ(wantPurelyInternal(global_constants::VF_HIGH + DBL_EPSILON), true);
  EXPECT_EQ(wantPurelyInternal(1.0), true);
}

TEST(Helper, wantPurelyExternal) {
  EXPECT_EQ(wantPurelyExternal(0.0), true);
  EXPECT_EQ(wantPurelyExternal(global_constants::VF_LOW - DBL_EPSILON), true);
  EXPECT_EQ(wantPurelyExternal(global_constants::VF_LOW + DBL_EPSILON), false);
  EXPECT_EQ(wantPurelyExternal(1.0), false);
  EXPECT_EQ(wantPurelyExternal(0.8), false);
}

TEST(Helper, sort3Ascending) {
  double items[3] = {68.2, -4.0, 3.8};
  sort3Ascending(items);
  EXPECT_DOUBLE_EQ(items[0], -4.0);
  EXPECT_DOUBLE_EQ(items[1], 3.8);
  EXPECT_DOUBLE_EQ(items[2], 68.2);
}

TEST(Helper, sort3Descending) {
  double items[3] = {68.2, -4.0, 3.8};
  sort3Descending(items);
  EXPECT_DOUBLE_EQ(items[0], 68.2);
  EXPECT_DOUBLE_EQ(items[1], 3.8);
  EXPECT_DOUBLE_EQ(items[2], -4.0);
}

TEST(Helper, sort3AscendingTracked) {
  auto items = std::array<double, 3>{{68.2, -4.0, 3.8}};
  auto item_loc = std::array<int, 3>{{0, 1, 2}};
  sort3AscendingTracked(items.data(), item_loc.data());
  EXPECT_DOUBLE_EQ(items[0], -4.0);
  EXPECT_DOUBLE_EQ(items[1], 3.8);
  EXPECT_DOUBLE_EQ(items[2], 68.2);
  EXPECT_EQ(item_loc[0], 1);
  EXPECT_EQ(item_loc[1], 2);
  EXPECT_EQ(item_loc[2], 0);
}

TEST(Helper, sort3DescendingTracked) {
  auto items = std::array<double, 3>{{68.2, -4.0, 3.8}};
  auto item_loc = std::array<int, 3>{{0, 1, 2}};
  sort3DescendingTracked(items.data(), item_loc.data());
  EXPECT_DOUBLE_EQ(items[0], 68.2);
  EXPECT_DOUBLE_EQ(items[1], 3.8);
  EXPECT_DOUBLE_EQ(items[2], -4.0);
  EXPECT_EQ(item_loc[0], 0);
  EXPECT_EQ(item_loc[1], 2);
  EXPECT_EQ(item_loc[2], 1);
}

TEST(Helper, sortAscendingBasedOnOtherArray) {
  std::array<double, 10> carried = {-100.0, 4.0, 2.0,  8.0,  6.0,
                                    2.0,    1.0, -8.0, -2.0, 5.0};
  std::array<int, 10> dictating = {
      62, 0, -2, 4, 8, 15, 6, 22, 1, 3,
  };
  sortAscendingBasedOnOtherArray(&carried, &dictating);
  std::array<double, 10> carried_sorted = {2.0, 4.0, -2.0, 5.0,  8.0,
                                           1.0, 6.0, 2.0,  -8.0, -100.0};
  std::array<int, 10> dictating_sorted = {-2, 0, 1, 3, 4, 6, 8, 15, 22, 62};
  for (UnsignedIndex_t n = 0; n < 10; ++n) {
    EXPECT_DOUBLE_EQ(carried[n], carried_sorted[n]);
    EXPECT_EQ(dictating[n], dictating_sorted[n]);
  }
}

TEST(Helper, sortDescendingBasedOnOtherArray) {
  std::array<double, 10> carried = {-100.0, 4.0, 2.0,  8.0,  6.0,
                                    2.0,    1.0, -8.0, -2.0, 5.0};
  std::array<int, 10> dictating = {
      62, 0, -2, 4, 8, 15, 6, 22, 1, 3,
  };
  sortDescendingBasedOnOtherArray(&carried, &dictating);
  std::array<double, 10> carried_sorted = {-100.0, -8.0, 2.0,  6.0, 1.0,
                                           8.0,    5.0,  -2.0, 4.0, 2.0};
  std::array<int, 10> dictating_sorted = {62, 22, 15, 8, 6, 4, 3, 1, 0, -2};

  for (UnsignedIndex_t n = 0; n < 10; ++n) {
    EXPECT_DOUBLE_EQ(carried[n], carried_sorted[n]);
    EXPECT_EQ(dictating[n], dictating_sorted[n]);
  }
}

}  // namespace
