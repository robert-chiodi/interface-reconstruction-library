// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/data_structures/lister.h"

#include "gtest/gtest.h"

namespace {

using namespace IRL;

TEST(ObjectCollections, Lister) {
  Lister<int> list_0;
  list_0.resize(4);
  list_0[0] = -4;
  list_0[1] = 2;
  list_0[2] = 68;
  list_0[3] = 12;
  list_0 += 3;
  list_0 += -50;
  int correct_values[6] = {-4, 2, 68, 12, 3, -50};
  EXPECT_EQ(list_0.size(), 6);
  int n = 0;
  for (const auto& number : list_0) {
    EXPECT_EQ(number, correct_values[n]);
    ++n;
  }

  Lister<int> list_1;
  list_1 += 15;
  list_1 += 6;
  list_0 += list_1;
  int correct_values_longer[8] = {-4, 2, 68, 12, 3, -50, 15, 6};
  EXPECT_EQ(list_0.size(), 8);
  n = 0;
  for (const auto& number : list_0) {
    EXPECT_EQ(number, correct_values_longer[n]);
    ++n;
  }
}

}  // namespace
