// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_BLOCK_OBJECT_ALLOCATION_H_
#define SRC_DATA_STRUCTURES_BLOCK_OBJECT_ALLOCATION_H_

#include "src/parameters/defined_types.h"

namespace IRL {

template <class ObjectType>
class BlockObjectAllocation {
 public:
  BlockObjectAllocation(void) = delete;

  explicit BlockObjectAllocation(const LargeOffsetIndex_t a_number_to_allocate);

  ObjectType* getObjectAtIndex(const LargeOffsetIndex_t a_object_index);

  explicit BlockObjectAllocation(const BlockObjectAllocation& other);

  explicit BlockObjectAllocation(BlockObjectAllocation&& other) noexcept;

  BlockObjectAllocation& operator=(const BlockObjectAllocation& other);

  BlockObjectAllocation& operator=(BlockObjectAllocation&& other) noexcept;

  ~BlockObjectAllocation(void);

 private:
  ObjectType* object_allocation_m;
  LargeOffsetIndex_t allocated_size_m;
};

}  // namespace IRL

#include "src/data_structures/block_object_allocation.tpp"

#endif  // SRC_DATA_STRUCTURES_BLOCK_OBJECT_ALLOCATION_H_
