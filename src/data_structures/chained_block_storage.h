// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_CHAINED_BLOCK_STORAGE_H_
#define SRC_DATA_STRUCTURES_CHAINED_BLOCK_STORAGE_H_

#include <algorithm>
#include <cassert>
#include <cstring>
#include <utility>
#include <tuple>
#include <vector>

#include <iostream>

#include "src/data_structures/iterator_through_bracket_operator.h"
#include "src/parameters/defined_types.h"

// Could add a free list and make this closer to a Slot Map data structure.

namespace IRL {

template <class ContainedType, UnsignedIndex_t kBlockSize>
class ChainedBlockStorage {
 public:
  using iterator = IteratorThroughBracketOperator<
      ChainedBlockStorage<ContainedType, kBlockSize>>;
  using const_iterator = ConstIteratorThroughBracketOperator<
      ChainedBlockStorage<ContainedType, kBlockSize>>;
  using value_t = ContainedType;

  ChainedBlockStorage(void);

  explicit ChainedBlockStorage(const UnsignedIndex_t a_initial_block_number);

  ContainedType& operator[](const UnsignedIndex_t a_index);

  const ContainedType& operator[](const UnsignedIndex_t a_index) const;

  ContainedType& getNextElement(void);
  ContainedType& getNextElement(const ContainedType& a_object);
  ContainedType& getNextElement(ContainedType&& a_object);	

  bool empty(void) const;

  UnsignedIndex_t size(void) const;

  // Note, calling resize will ruin the validity of any pointers held to the
  // objects. It will assume that the objects from 0 to size-1 will be
  // referenced/taken from ChainedBlockStorage::operator[]()
  void resize(const UnsignedIndex_t a_size);

  UnsignedIndex_t currentSupportedSize(void) const;

  // Will reserve atleast enough space for a_size,
  // while rounding up to the creation of a kBlockSize block.
  void reserve(const UnsignedIndex_t a_size);

  void freeObject(ContainedType* a_object_to_free);
  void deallocateMemory(void);

  // Copy constructor
  ChainedBlockStorage(const ChainedBlockStorage& a_rhs) noexcept;

  // Copy Assignment
  ChainedBlockStorage& operator=(const ChainedBlockStorage& a_rhs) noexcept;

  iterator begin(void);
  iterator end(void);
  const_iterator begin(void) const;
  const_iterator end(void) const;
  const_iterator cbegin(void) const;
  const_iterator cend(void) const;

 private:
  ContainedType* allocateBlock(void);
  void increaseToSize(const UnsignedIndex_t a_size);
  ContainedType& getNextFreeObject(void);
  void incrementNextFreeObject(void);
  bool noMoreFreeObjects(void) const;
  void addBlockToFreeList(const UnsignedIndex_t a_block_index);

  std::vector<ContainedType*> data_blocks_m;

  struct FreeObjectLocation{

	  FreeObjectLocation(const std::size_t a_block, const std::size_t a_index)
	  : block(a_block), index(a_index) {}

	  std::size_t block;
	  std::size_t index;
  };
  FreeObjectLocation next_free_object_m;
};

}  // namespace IRL
#include "src/data_structures/chained_block_storage.tpp"

#endif  // SRC_DATA_STRUCTURES_CHAINED_BLOCK_STORAGE_H_
