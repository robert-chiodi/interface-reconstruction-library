// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_DATA_STRUCTURES_CHAINED_BLOCK_STORAGE_H_
#define IRL_DATA_STRUCTURES_CHAINED_BLOCK_STORAGE_H_

#include <algorithm>
#include <cassert>
#include <cstring>
#include <tuple>
#include <utility>
#include <vector>

#include <iostream>

#include "irl/data_structures/iterator_through_bracket_operator.h"
#include "irl/parameters/defined_types.h"

// Could add a free list and make this closer to a Slot Map data structure.

namespace IRL {

template <class ContainedType, UnsignedIndex_t kBlockSize>
class ChainedBlockStorage {
public:
  using value_t = ContainedType;

  ChainedBlockStorage(void);

  ContainedType &operator[](const UnsignedIndex_t a_index);

  const ContainedType &operator[](const UnsignedIndex_t a_index) const;

  ContainedType &getNextElement(void);
  ContainedType &getNextElement(const ContainedType &a_object);
  ContainedType &getNextElement(ContainedType &&a_object);

  UnsignedIndex_t size(void) const;

  // Note, calling resize will ruin the validity of any pointers held to the
  // objects. It will assume that the objects from 0 to size-1 will be
  // referenced/taken from ChainedBlockStorage::operator[]()
  void resize(const UnsignedIndex_t a_size);

  UnsignedIndex_t currentSupportedSize(void) const;

  // Will reserve atleast enough space for a_size,
  // while rounding up to the creation of a kBlockSize block.
  void reserve(const UnsignedIndex_t a_size);

  void deallocateMemory(void);

  // Copy constructor
  ChainedBlockStorage(const ChainedBlockStorage &a_rhs) noexcept;

  // Copy Assignment
  ChainedBlockStorage &operator=(const ChainedBlockStorage &a_rhs) noexcept;

  // Destructor
  ~ChainedBlockStorage(void);

private:
  ContainedType *allocateBlock(void);
  void increaseToSize(const UnsignedIndex_t a_size);
  ContainedType &getNextFreeObject(void);
  void incrementNextFreeObject(void);
  bool noMoreFreeObjects(void) const;

  std::vector<ContainedType *> data_blocks_m;

  ContainedType *next_free_object_m;
  UnsignedIndex_t open_block_m;
};

} // namespace IRL
#include "irl/data_structures/chained_block_storage.tpp"

#endif // SRC_DATA_STRUCTURES_CHAINED_BLOCK_STORAGE_H_
