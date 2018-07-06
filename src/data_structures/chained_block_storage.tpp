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
    : data_blocks_m(), next_free_object_m(0, 0) {}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ChainedBlockStorage<ContainedType, kBlockSize>::ChainedBlockStorage(
    const UnsignedIndex_t a_initial_block_number)
    : data_blocks_m(a_initial_block_number), next_free_object_m(0, 0) {
  for (UnsignedIndex_t n = 0; n < a_initial_block_number; ++n) {
    data_blocks_m[n] = this->allocateBlock();
  }
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ContainedType& ChainedBlockStorage<ContainedType, kBlockSize>::operator[](
    const UnsignedIndex_t a_index) {
  assert(a_index < this->currentSupportedSize());
  const auto block_index = a_index / kBlockSize;
  return data_blocks_m[block_index][a_index - block_index * kBlockSize];
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
const ContainedType& ChainedBlockStorage<ContainedType, kBlockSize>::operator[](
    const UnsignedIndex_t a_index) const {
  assert(a_index < this->currentSupportedSize());
  const auto block_index = a_index / kBlockSize;
  return data_blocks_m[block_index][a_index - block_index * kBlockSize];
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ContainedType&
ChainedBlockStorage<ContainedType, kBlockSize>::getNextFreeObject(void) {
  assert(next_free_object_m.block < data_blocks_m.size());
  assert(next_free_object_m.index < kBlockSize);
  ContainedType& object_to_return =
      data_blocks_m[next_free_object_m.block][next_free_object_m.index];
  this->incrementNextFreeObject();
  return object_to_return;
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
void ChainedBlockStorage<ContainedType, kBlockSize>::incrementNextFreeObject(
    void) {
  ++next_free_object_m.index;
  if (next_free_object_m.index == kBlockSize) {
    ++next_free_object_m.block;
    assert(next_free_object_m.block <= data_blocks_m.size());
    next_free_object_m.index = 0;
  }
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
bool ChainedBlockStorage<ContainedType, kBlockSize>::noMoreFreeObjects(
    void) const {
  return next_free_object_m.block == data_blocks_m.size();
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ContainedType& ChainedBlockStorage<ContainedType, kBlockSize>::getNextElement(
    void) {
  if (this->noMoreFreeObjects()) {
    this->increaseToSize(this->size() + 1);
  }
  ContainedType& object_to_return = this->getNextFreeObject();
  object_to_return = ContainedType();
  return object_to_return;
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ContainedType& ChainedBlockStorage<ContainedType, kBlockSize>::getNextElement(
    const ContainedType& a_object) {
  if (this->noMoreFreeObjects()) {
    this->increaseToSize(this->size() + 1);
  }
  ContainedType& object_to_return = this->getNextFreeObject();
  object_to_return = a_object;
  return object_to_return;
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ContainedType& ChainedBlockStorage<ContainedType, kBlockSize>::getNextElement(
    ContainedType&& a_object) {
  if (this->noMoreFreeObjects()) {
    this->increaseToSize(this->size() + 1);
  }
  ContainedType& object_to_return = this->getNextFreeObject();
  object_to_return = std::move(a_object);
  return object_to_return;
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
bool ChainedBlockStorage<ContainedType, kBlockSize>::empty(void) const {
  return this->size() == 0;
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
UnsignedIndex_t ChainedBlockStorage<ContainedType, kBlockSize>::size(
    void) const {
  return static_cast<UnsignedIndex_t>(next_free_object_m.block) * kBlockSize +
         static_cast<UnsignedIndex_t>(next_free_object_m.index);
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
void ChainedBlockStorage<ContainedType, kBlockSize>::resize(
    const UnsignedIndex_t a_size) {
  // Increase to sufficient size
  this->increaseToSize(a_size);
  next_free_object_m.block = a_size / kBlockSize;
  next_free_object_m.index = a_size % kBlockSize;
  assert(this->size() == a_size);
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
UnsignedIndex_t
ChainedBlockStorage<ContainedType, kBlockSize>::currentSupportedSize(
    void) const {
  return static_cast<UnsignedIndex_t>(data_blocks_m.size()) * kBlockSize;
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
void ChainedBlockStorage<ContainedType, kBlockSize>::reserve(
    const UnsignedIndex_t a_size) {
  this->increaseToSize(a_size);
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
void ChainedBlockStorage<ContainedType, kBlockSize>::freeObject(
    ContainedType* a_object_to_free) {
  std::cout << "Freeing of objects is currently disabled" << std::endl;
  std::exit(-1);
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
void ChainedBlockStorage<ContainedType, kBlockSize>::deallocateMemory(void) {
  for (auto& elem : data_blocks_m) {
    delete[] elem;
  }
  data_blocks_m.clear();
  next_free_object_m = FreeObjectLocation(0, 0);
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ChainedBlockStorage<ContainedType, kBlockSize>::ChainedBlockStorage(
    const ChainedBlockStorage& a_rhs) noexcept
    : data_blocks_m(), next_free_object_m(0, 0) {
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
ChainedBlockStorage<ContainedType, kBlockSize>&
ChainedBlockStorage<ContainedType, kBlockSize>::operator=(
    const ChainedBlockStorage& a_rhs) noexcept {
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
typename ChainedBlockStorage<ContainedType, kBlockSize>::iterator
ChainedBlockStorage<ContainedType, kBlockSize>::begin(void) {
  return iterator(*this, 0);
}
template <class ContainedType, UnsignedIndex_t kBlockSize>
typename ChainedBlockStorage<ContainedType, kBlockSize>::iterator
ChainedBlockStorage<ContainedType, kBlockSize>::end(void) {
  return iterator(*this, this->size());
}
template <class ContainedType, UnsignedIndex_t kBlockSize>
typename ChainedBlockStorage<ContainedType, kBlockSize>::const_iterator
ChainedBlockStorage<ContainedType, kBlockSize>::begin(void) const {
  return this->cbegin();
}
template <class ContainedType, UnsignedIndex_t kBlockSize>
typename ChainedBlockStorage<ContainedType, kBlockSize>::const_iterator
ChainedBlockStorage<ContainedType, kBlockSize>::end(void) const {
  return this->cend();
}
template <class ContainedType, UnsignedIndex_t kBlockSize>
typename ChainedBlockStorage<ContainedType, kBlockSize>::const_iterator
ChainedBlockStorage<ContainedType, kBlockSize>::cbegin(void) const {
  return const_iterator(*this, 0);
}
template <class ContainedType, UnsignedIndex_t kBlockSize>
typename ChainedBlockStorage<ContainedType, kBlockSize>::const_iterator
ChainedBlockStorage<ContainedType, kBlockSize>::cend(void) const {
  return const_iterator(*this, this->size());
}

template <class ContainedType, UnsignedIndex_t kBlockSize>
ContainedType* ChainedBlockStorage<ContainedType, kBlockSize>::allocateBlock(
    void) {
  return static_cast<ContainedType*>(::operator new(
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
  if (new_blocks_needed > 0 && this->noMoreFreeObjects()) {
    next_free_object_m.block = data_blocks_m.size();
    next_free_object_m.index = 0;
  }
  for (int n = 0; n < new_blocks_needed; ++n) {
    data_blocks_m.push_back(this->allocateBlock());
  }
}

}  // namespace IRL

#endif  // SRC_DATA_STRUCTURES_CHAINED_BLOCK_STORAGE_TPP_
