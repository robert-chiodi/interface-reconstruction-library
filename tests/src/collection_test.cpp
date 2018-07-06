// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/data_structures/collection.h"

#include "gtest/gtest.h"

namespace {

using namespace IRL;

TEST(ObjectCollections, Collection) {
  Collection<int> collection;
  collection.resize(4);
  collection[0] = -4;
  collection[1] = 2;
  collection[2] = 68;
  collection[3] = 12;
  collection.push_back(17);
  collection.emplace_back(-8);
  int correct_values[6] = {-4, 2, 68, 12, 17, -8};
  EXPECT_EQ(collection.size(), 6);
  int n = 0;
  for (const auto& number : collection) {
    EXPECT_EQ(number, correct_values[n]);
    ++n;
  }
  collection[1]++;
  correct_values[1]++;
  n = 0;
  for (const auto& number : collection) {
    EXPECT_EQ(number, correct_values[n]);
    ++n;
  }
}
}  // namespace
