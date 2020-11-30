// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/generic_cutting/general/encountered_pair_list.h"

#include "gtest/gtest.h"

namespace {

using namespace IRL;

TEST(EncounteredPairList, EncounteredPairList) {
  EncounteredPairList list;

  std::array<UnsignedIndex_t, 6> first_row{{16, 14, 14, 48, 52, 1}};
  std::array<UnsignedIndex_t, 6> second_row{{22, 7, 6, 920, 16, 8220}};
  for (UnsignedIndex_t n = 0; n < first_row.size(); ++n) {
    list.addPair(first_row[n], second_row[n]);
  }

  list.printOut();

  for (UnsignedIndex_t n = 0; n < first_row.size(); ++n) {
    EXPECT_TRUE(list.isPairPresent(first_row[n], second_row[n]))
        << "Pair missing is " << first_row[n] << " " << second_row[n]
        << std::endl;
    EXPECT_TRUE(list.isPairPresent(second_row[n], first_row[n]))
        << "Pair missing is " << first_row[n] << " " << second_row[n]
        << std::endl;
  }
  EXPECT_FALSE(list.isPairPresent(86, 87));
  EXPECT_FALSE(list.isPairPresent(14, 23));

  auto initial_state = list.saveState();
  std::array<UnsignedIndex_t, 3> add_first_row{{86, 14, 6}};
  std::array<UnsignedIndex_t, 3> add_second_row{{87, 23, 22}};
  for (UnsignedIndex_t n = 0; n < add_first_row.size(); ++n) {
    list.addPair(add_first_row[n], add_second_row[n]);
  }
  std::cout << "List with additional members" << std::endl;
  list.printOut();
  for (UnsignedIndex_t n = 0; n < first_row.size(); ++n) {
    EXPECT_TRUE(list.isPairPresent(first_row[n], second_row[n]))
        << "Pair missing is " << first_row[n] << " " << second_row[n]
        << std::endl;
    EXPECT_TRUE(list.isPairPresent(second_row[n], first_row[n]))
        << "Pair missing is " << first_row[n] << " " << second_row[n]
        << std::endl;
  }
  for (UnsignedIndex_t n = 0; n < add_first_row.size(); ++n) {
    EXPECT_TRUE(list.isPairPresent(add_first_row[n], add_second_row[n]))
        << "Pair missing is " << add_first_row[n] << " " << add_second_row[n]
        << std::endl;
    EXPECT_TRUE(list.isPairPresent(add_second_row[n], add_first_row[n]))
        << "Pair missing is " << add_first_row[n] << " " << add_second_row[n]
        << std::endl;
  }
  list.resetToState(initial_state);
  for (UnsignedIndex_t n = 0; n < first_row.size(); ++n) {
    EXPECT_TRUE(list.isPairPresent(first_row[n], second_row[n]))
        << "Pair missing is " << first_row[n] << " " << second_row[n]
        << std::endl;
    EXPECT_TRUE(list.isPairPresent(second_row[n], first_row[n]))
        << "Pair missing is " << first_row[n] << " " << second_row[n]
        << std::endl;
  }
  for (UnsignedIndex_t n = 0; n < add_first_row.size(); ++n) {
    EXPECT_FALSE(list.isPairPresent(add_first_row[n], add_second_row[n]))
        << "Pair should not exist " << add_first_row[n] << " "
        << add_second_row[n] << std::endl;
    EXPECT_FALSE(list.isPairPresent(add_second_row[n], add_first_row[n]))
        << "Pair should not exist " << add_first_row[n] << " "
        << add_second_row[n] << std::endl;
  }
}

} // namespace
