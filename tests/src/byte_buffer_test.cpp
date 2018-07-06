// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/helpers/byte_buffer.h"

#include "gtest/gtest.h"

#include <array>

#include "src/parameters/defined_types.h"

namespace {

using namespace IRL;

TEST(Serializer, ByteBuffer) {
  // Set known int and double arrays
  std::array<int, 2> int_array = {14, -4};
  std::array<double, 6> double_array = {-4.8, 6.2, -2.9, 100.43, 42.96, -18.88};
  // Set single int/double to tack onto end of buffer
  int single_int = 22;
  double single_double = 96.089;

  // Construct and pack buffer
  ByteBuffer buffer;
  buffer.pack(int_array.data(), int_array.size());
  buffer.pack(double_array.data(), double_array.size());
  buffer.pack(&single_int, 1);
  buffer.pack(&single_double, 1);

  // Reset buffer pointer to read from beginning
  buffer.resetBufferPointer();

  // Allocate receiving arrays and unpack into them
  std::array<int, 3> recv_int;
  std::array<double, 7> recv_double;
  buffer.unpack(recv_int.data(), 2);
  buffer.unpack(recv_double.data(), 6);
  buffer.unpack(&recv_int[2], 1);
  buffer.unpack(&recv_double[6], 1);

  // Check unpacked values are same as packed
  for (UnsignedIndex_t n = 0; n < 2; ++n) {
    EXPECT_DOUBLE_EQ(recv_int[n], int_array[n]);
  }
  for (UnsignedIndex_t n = 0; n < 6; ++n) {
    EXPECT_DOUBLE_EQ(recv_double[n], double_array[n]);
  }
  EXPECT_DOUBLE_EQ(recv_int[2], single_int);
  EXPECT_DOUBLE_EQ(recv_double[6], single_double);
}
}  // namespace
