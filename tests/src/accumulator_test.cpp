// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/data_structures/accumulator.h"

#include "gtest/gtest.h"

namespace {

using namespace IRL;

TEST(ObjectCollections, Accumulator) {
  Accumulator<int> accumulator_0;
  accumulator_0.resize(4);
  accumulator_0[0] = -4;
  accumulator_0[1] = 2;
  accumulator_0[2] = 68;
  accumulator_0[3] = 12;
  accumulator_0.push_back(17);
  accumulator_0.emplace_back(-8);
  accumulator_0[6] = 14;
  int correct_values[7] = {-4, 2, 68, 12, 17, -8, 14};
  EXPECT_EQ(accumulator_0.size(), 7);
  int n = 0;
  for (const auto& number : accumulator_0) {
    EXPECT_EQ(number, correct_values[n]);
    ++n;
  }

  Accumulator<int> accumulator_1;
  accumulator_1.resize(4);
  accumulator_1[0] = -4;
  accumulator_1[1] = 2;
  accumulator_1[2] = 68;
  accumulator_1[3] = 12;
  accumulator_1.push_back(17);
  accumulator_1.emplace_back(-8);
  accumulator_1[6] = 14;
  accumulator_1[7] = 3;
  accumulator_1[8] = 8;
  accumulator_1[9] = 6;
  int correct_values_longer[10] = {-8, 4, 136, 24, 34, -16, 28, 3, 8, 6};
  accumulator_0 += accumulator_1;
  EXPECT_EQ(accumulator_0.size(), 10);
  n = 0;
  for (const auto& number : accumulator_0) {
    EXPECT_EQ(number, correct_values_longer[n]);
    ++n;
  }
  Accumulator<int> accumulator_2;
  accumulator_2[1] = 1;
  accumulator_2[2] = 1;
  accumulator_2[0] = 1;
  accumulator_2[3] = 1;
  correct_values_longer[0]++;
  correct_values_longer[1]++;
  correct_values_longer[2]++;
  correct_values_longer[3]++;
  accumulator_0 += accumulator_2;
  EXPECT_EQ(accumulator_0.size(), 10);
  n = 0;
  for (const auto& number : accumulator_0) {
    EXPECT_EQ(number, correct_values_longer[n]);
    ++n;
  }
}
}  // namespace
