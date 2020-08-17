// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/helpers/serializer.h"

#include "gtest/gtest.h"

#include <array>
#include <limits>
#include <random>

#include "src/helpers/byte_buffer.h"
#include "src/parameters/defined_types.h"

namespace {

namespace {

std::uniform_real_distribution<double> getDistribution(double a_dummy) {
  return std::uniform_real_distribution<double>(-100000.0, 100000.0);
}

std::uniform_int_distribution<int> getDistribution(int a_dummy) {
  return std::uniform_int_distribution<int>(-100000, 100000);
}

}  // namespace

template <class ElementType, IRL::UnsignedIndex_t kAmount>
class RandomNumbers {
 public:
  RandomNumbers(void) {
    std::random_device
        rd;  // Get a random seed from the OS entropy device, or whatever
    std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                                // generator and seed it with entropy.
    ElementType i = 0;
    auto random_number = getDistribution(i);
    for (auto& elem : numbers_m) {
      elem = random_number(eng);
    }
  }

  ElementType& operator[](const IRL::UnsignedIndex_t a_index) {
    assert(a_index < kAmount);
    return numbers_m[a_index];
  }

  ElementType* data(void) { return numbers_m.data(); }
  const ElementType* data(void) const { return numbers_m.data(); }

  void serialize(IRL::ByteBuffer* buffer_to_pack) const {
    buffer_to_pack->pack(this->data(), numbers_m.size());
  }

  void unpackSerialized(IRL::ByteBuffer* buffer_to_unpack) {
    buffer_to_unpack->unpack(this->data(), numbers_m.size());
  }

 private:
  std::array<ElementType, kAmount> numbers_m;
};

using namespace IRL;

TEST(Serializer, ByteBufferRandom) {
  static constexpr int amount_of_numbers = 100;
  // Allocate arrays of random doubles and int's
  RandomNumbers<double, amount_of_numbers> random_doubles_0;
  RandomNumbers<int, amount_of_numbers> random_ints_0;
  RandomNumbers<double, amount_of_numbers> random_doubles_1;
  RandomNumbers<int, amount_of_numbers> random_ints_1;

  // Construct ByteBuffer
  ByteBuffer buffer;

  // Use RandomNumbers::serialize to pack this buffer
  serializeAndPack(random_doubles_0, &buffer);
  serializeAndPack(random_ints_0, &buffer);
  serializeAndPack(random_doubles_1, &buffer);
  serializeAndPack(random_ints_1, &buffer);

  // Reset ByteBuffer pointer to start for reading
  buffer.resetBufferPointer();

  // Create new RandomNumbers array and unpack into them.
  RandomNumbers<double, amount_of_numbers> recv_random_doubles_0;
  RandomNumbers<int, amount_of_numbers> recv_random_ints_0;
  RandomNumbers<double, amount_of_numbers> recv_random_doubles_1;
  RandomNumbers<int, amount_of_numbers> recv_random_ints_1;

  unpackAndStore(&recv_random_doubles_0, &buffer);
  unpackAndStore(&recv_random_ints_0, &buffer);
  unpackAndStore(&recv_random_doubles_1, &buffer);
  unpackAndStore(&recv_random_ints_1, &buffer);

  // Check packed and unpacked values are same.
  for (UnsignedIndex_t i = 0; i < amount_of_numbers; ++i) {
    EXPECT_DOUBLE_EQ(recv_random_doubles_0[i], random_doubles_0[i]);
    EXPECT_DOUBLE_EQ(recv_random_doubles_1[i], random_doubles_1[i]);
  }
  for (UnsignedIndex_t i = 0; i < amount_of_numbers; ++i) {
    EXPECT_EQ(recv_random_ints_0[i], random_ints_0[i]);
    EXPECT_EQ(recv_random_ints_1[i], random_ints_1[i]);
  }
}

}  // namespace
