// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_STACK_VECTOR_H_
#define SRC_DATA_STRUCTURES_STACK_VECTOR_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>

#include "src/parameters/defined_types.h"

namespace IRL {

template <class ObjectType, UnsignedIndex_t kMaxSize>
class StackVector {
 public:
  using iterator = typename std::array<ObjectType, kMaxSize>::iterator;
  using const_iterator =
      typename std::array<ObjectType, kMaxSize>::const_iterator;

  StackVector(void);

  template <UnsignedIndex_t kOtherMaxSize>
  explicit StackVector(const StackVector<ObjectType, kOtherMaxSize>& other);

  template <UnsignedIndex_t kOtherMaxSize>
  StackVector& operator=(const StackVector<ObjectType, kOtherMaxSize>& other);

  ObjectType& operator[](const UnsignedIndex_t a_index);

  const ObjectType& operator[](const UnsignedIndex_t a_index) const;

  UnsignedIndex_t size(void) const;

  UnsignedIndex_t maxSize(void) const;

  void resize(const UnsignedIndex_t a_count);
  void assign(const UnsignedIndex_t a_count, const ObjectType& a_value);

  ObjectType& front(void);
  const ObjectType& front(void) const;

  ObjectType& back(void);

  const ObjectType& back(void) const;

  bool empty(void) const;
  bool full(void) const;

  void pop_back(void);

  void push_back(const ObjectType& a_object);

  void erase(iterator a_pos);

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  ~StackVector(void) = default;

 private:
  void incrementSize(void);

  void decrementSize(void);

  std::array<ObjectType, kMaxSize> object_container_m;
  UnsignedIndex_t size_m;
};

template <class ObjectType, UnsignedIndex_t kMaxSize>
inline std::ostream& operator<<(
    std::ostream& out, const StackVector<ObjectType, kMaxSize>& a_stack_vector);

}  // namespace IRL

#include "src/data_structures/stack_vector.tpp"

#endif  // SRC_DATA_STRUCTURES_STACK_VECTOR_H_
