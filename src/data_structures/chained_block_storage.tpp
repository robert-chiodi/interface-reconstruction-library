// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_CHAINED_BLOCK_STORAGE_TPP_
#define SRC_DATA_STRUCTURES_CHAINED_BLOCK_STORAGE_TPP_

#include <iostream>

namespace IRL {

template <class ContainedType, UnsignedIndex_t kBlockSize>
ChainedBlockStorage<ContainedType, kBlockSize>::ChainedBlockStorage(void)
    : data_blocks_m(), next_free_object_m(nullptr),
      open_block_m(static_cast<UnsignedIndex_t>(-1)) {}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ContainedType &ChainedBlockStorage<ContainedType, kBlockSize>::
operator[](const UnsignedIndex_t a_index) {
  assert(a_index < this->size());
  const auto block_index = a_index / kBlockSize;
  return data_blocks_m[block_index][a_index - block_index * kBlockSize];
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
const ContainedType &ChainedBlockStorage<ContainedType, kBlockSize>::
operator[](const UnsignedIndex_t a_index) const {
  assert(a_index < this->size());
  const auto block_index = a_index / kBlockSize;
  return data_blocks_m[block_index][a_index - block_index * kBlockSize];
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ContainedType &
ChainedBlockStorage<ContainedType, kBlockSize>::getNextFreeObject(void) {
  if (next_free_object_m - data_blocks_m[open_block_m] == kBlockSize) {
    ++open_block_m;
    if (open_block_m == data_blocks_m.size()) {
      data_blocks_m.push_back(allocateBlock());
    }
    next_free_object_m = data_blocks_m[open_block_m];
  }
  return *(next_free_object_m++);
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ContainedType &
ChainedBlockStorage<ContainedType, kBlockSize>::getNextElement(void) {
  ContainedType &object_to_return = this->getNextFreeObject();
  object_to_return = ContainedType();
  return object_to_return;
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ContainedType &ChainedBlockStorage<ContainedType, kBlockSize>::getNextElement(
    const ContainedType &a_object) {
  ContainedType &object_to_return = this->getNextFreeObject();
  object_to_return = a_object;
  return object_to_return;
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ContainedType &ChainedBlockStorage<ContainedType, kBlockSize>::getNextElement(
    ContainedType &&a_object) {
  ContainedType &object_to_return = this->getNextFreeObject();
  object_to_return = std::move(a_object);
  return object_to_return;
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
UnsignedIndex_t
ChainedBlockStorage<ContainedType, kBlockSize>::size(void) const {
  if (data_blocks_m.empty()) {
    return 0;
  } else {
    return static_cast<UnsignedIndex_t>(open_block_m) * kBlockSize +
           static_cast<UnsignedIndex_t>(next_free_object_m -
                                        data_blocks_m[open_block_m]);
  }
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
void ChainedBlockStorage<ContainedType, kBlockSize>::resize(
    const UnsignedIndex_t a_size) {
  // Increase to sufficient size
  this->increaseToSize(a_size);
  open_block_m = a_size / kBlockSize;
  next_free_object_m = data_blocks_m[open_block_m] + a_size % kBlockSize;
  assert(this->size() == a_size);
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
UnsignedIndex_t
ChainedBlockStorage<ContainedType, kBlockSize>::currentSupportedSize(
    void) const {
  return static_cast<UnsignedIndex_t>(data_blocks_m.size()) * kBlockSize;
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
void ChainedBlockStorage<ContainedType, kBlockSize>::deallocateMemory(void) {
  for (auto &elem : data_blocks_m) {
    delete[] elem;
  }
  data_blocks_m.clear();
  next_free_object_m = nullptr;
  open_block_m = static_cast<UnsignedIndex_t>(-1);
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ChainedBlockStorage<ContainedType, kBlockSize>::ChainedBlockStorage(
    const ChainedBlockStorage &a_rhs) noexcept
    : data_blocks_m() {
  this->resize(a_rhs.size());
  assert(this->currentSupportedSize() >= a_rhs.size());
  std::size_t amount_to_copy =
      static_cast<std::size_t>(a_rhs.size()) * sizeof(ContainedType);
  static constexpr std::size_t block_size =
      static_cast<std::size_t>(kBlockSize) * sizeof(ContainedType);
  const std::size_t blocks_to_copy = amount_to_copy / block_size;
  for (std::size_t n = 0; n < blocks_to_copy; ++n) {
    std::memcpy(data_blocks_m[n], a_rhs.data_blocks_m[n], block_size);
  }
  assert(static_cast<int>(amount_to_copy) -
                 static_cast<int>(blocks_to_copy * block_size) >=
             0 ||
         (data_blocks_m.size() > blocks_to_copy &&
          a_rhs.data_blocks_m.size() > blocks_to_copy));
  std::memcpy(data_blocks_m[blocks_to_copy],
              a_rhs.data_blocks_m[blocks_to_copy],
              amount_to_copy - blocks_to_copy * block_size);
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ChainedBlockStorage<ContainedType, kBlockSize> &
ChainedBlockStorage<ContainedType, kBlockSize>::
operator=(const ChainedBlockStorage &a_rhs) noexcept {
  if (this != &a_rhs) {
    this->resize(a_rhs.size());
    assert(this->currentSupportedSize() >= a_rhs.size());
    std::size_t amount_to_copy =
        static_cast<std::size_t>(a_rhs.size()) * sizeof(ContainedType);
    static constexpr std::size_t block_size =
        static_cast<std::size_t>(kBlockSize) * sizeof(ContainedType);
    const std::size_t blocks_to_copy = amount_to_copy / block_size;
    for (std::size_t n = 0; n < blocks_to_copy; ++n) {
      std::memcpy(data_blocks_m[n], a_rhs.data_blocks_m[n], block_size);
    }
    assert(static_cast<int>(amount_to_copy) -
                   static_cast<int>(blocks_to_copy * block_size) >=
               0 ||
           (data_blocks_m.size() > blocks_to_copy &&
            a_rhs.data_blocks_m.size() > blocks_to_copy));
    std::memcpy(data_blocks_m[blocks_to_copy],
                a_rhs.data_blocks_m[blocks_to_copy],
                amount_to_copy - blocks_to_copy * block_size);
  }
  return (*this);
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ContainedType *
ChainedBlockStorage<ContainedType, kBlockSize>::allocateBlock(void) {
  return static_cast<ContainedType *>(::operator new(
      static_cast<std::size_t>(kBlockSize) * sizeof(ContainedType)));
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
void ChainedBlockStorage<ContainedType, kBlockSize>::increaseToSize(
    const UnsignedIndex_t a_size) {
  const int remaining_space =
      static_cast<int>(this->currentSupportedSize()) - static_cast<int>(a_size);
  const int new_blocks_needed =
      remaining_space < 0 ? -remaining_space / static_cast<int>(kBlockSize) + 1
                          : 0;
  for (int n = 0; n < new_blocks_needed; ++n) {
    data_blocks_m.push_back(this->allocateBlock());
  }
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ChainedBlockStorage<ContainedType, kBlockSize>::~ChainedBlockStorage(void) {
  for (auto &block : data_blocks_m) {
    ::operator delete(block);
    block = nullptr;
  }
  open_block_m = static_cast<UnsignedIndex_t>(-1);
  next_free_object_m = nullptr;
  data_blocks_m.resize(0);
}

} // namespace IRL

#endif // SRC_DATA_STRUCTURES_CHAINED_BLOCK_STORAGE_TPP_
