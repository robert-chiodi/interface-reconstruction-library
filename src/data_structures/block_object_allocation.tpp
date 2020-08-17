// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_BLOCK_OBJECT_ALLOCATION_TPP_
#define SRC_DATA_STRUCTURES_BLOCK_OBJECT_ALLOCATION_TPP_

#include <cassert>
#include <cstring>
#include <limits>

namespace IRL {

template <class ObjectType>
BlockObjectAllocation<ObjectType>::BlockObjectAllocation(
    const LargeOffsetIndex_t a_number_to_allocate)
    : allocated_size_m{a_number_to_allocate} {
  object_allocation_m = new ObjectType[allocated_size_m];
}

template <class ObjectType>
ObjectType* BlockObjectAllocation<ObjectType>::getObjectAtIndex(
    const LargeOffsetIndex_t a_object_index) {
  assert(a_object_index < allocated_size_m);
  return &object_allocation_m[a_object_index];
}

template <class ObjectType>
BlockObjectAllocation<ObjectType>::BlockObjectAllocation(
    const BlockObjectAllocation& other) {
  this->~BlockObjectAllocation();
  allocated_size_m = other.allocated_size_m;
  object_allocation_m = new ObjectType[allocated_size_m];
  std::memcpy(object_allocation_m, other.object_allocation_m,
              static_cast<std::size_t>(allocated_size_m) * sizeof(ObjectType));
}

template <class ObjectType>
BlockObjectAllocation<ObjectType>::BlockObjectAllocation(
    BlockObjectAllocation&& other) noexcept {
  this->~BlockObjectAllocation();
  allocated_size_m = other.allocated_size_m;
  object_allocation_m = other.object_allocation_m;
  other.object_allocation_m = nullptr;
  other.allocated_size_m = 0;
}

template <class ObjectType>
BlockObjectAllocation<ObjectType>& BlockObjectAllocation<ObjectType>::operator=(
    const BlockObjectAllocation& other) {
  if (this != &other) {
    this->~BlockObjectAllocation();
    allocated_size_m = other.allocated_size_m;
    object_allocation_m = new ObjectType[allocated_size_m];
    std::memcpy(
        object_allocation_m, other.object_allocation_m,
        static_cast<std::size_t>(allocated_size_m) * sizeof(ObjectType));
  }
  return (*this);
}

template <class ObjectType>
BlockObjectAllocation<ObjectType>& BlockObjectAllocation<ObjectType>::operator=(
    BlockObjectAllocation&& other) noexcept {
  if (this != other) {
    this->~BlockObjectAllocation();
    allocated_size_m = other.allocated_size_m;
    object_allocation_m = other.object_allocation_m;
    other.object_allocation_m = nullptr;
    other.allocated_size_m = 0;
  }
  return (*this);
}

template <class ObjectType>
BlockObjectAllocation<ObjectType>::~BlockObjectAllocation(void) {
  delete[] object_allocation_m;
  object_allocation_m = nullptr;
  allocated_size_m = 0;
}

}  // namespace IRL

#endif  // SRC_DATA_STRUCTURES_BLOCK_OBJECT_ALLOCATION_TPP_
