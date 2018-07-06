// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/data_structures/stack_vector.h"

#include "gtest/gtest.h"

#include "src/parameters/defined_types.h"

namespace {

using namespace IRL;

TEST(StackVector, erase) {
  StackVector<UnsignedIndex_t, 10> test_vector;
  test_vector.resize(4);
  ASSERT_EQ(test_vector.size(), 4);
  test_vector[0] = 4;
  test_vector[1] = 8;
  test_vector[2] = 16;
  test_vector[3] = 222;
  test_vector.erase(test_vector.begin() + 1);
  EXPECT_EQ(test_vector.size(), 3);
  EXPECT_EQ(test_vector[0], 4);
  EXPECT_EQ(test_vector[1], 16);
  EXPECT_EQ(test_vector[2], 222);
}

}  // namespace
