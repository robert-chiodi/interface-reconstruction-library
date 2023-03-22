// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_DATA_STRUCTURES_SMALL_VECTOR_H_
#define IRL_DATA_STRUCTURES_SMALL_VECTOR_H_

#ifndef IRL_NO_ABSL
#include "absl/container/inlined_vector.h"
namespace IRL {
template <class T, std::size_t BuffSizeInElements>
using SmallVector =
    absl::InlinedVector<T, BuffSizeInElements, std::allocator<T>>;
}  // namespace IRL
#endif  // #ifndef IRL_NO_ABSL

#ifdef IRL_NO_ABSL
#include <cstddef>
#include <initializer_list>
#include <type_traits>
#include <utility>
#include <vector>

#include "irl/data_structures/short_alloc.h"
#include "irl/parameters/defined_types.h"

namespace IRL {
template <class T, std::size_t BuffSizeInElements>
class SmallVector {
  using BaseVector = std::vector<
      T, short_alloc::short_alloc<
             T, static_cast<std::size_t>(BuffSizeInElements) * sizeof(T),
             alignof(T)>>;

  static_assert(!std::is_same<T, bool>::value,
                "SmallVector does not support bool.");

 public:
  using iterator = typename BaseVector::iterator;
  using const_iterator = typename BaseVector::const_iterator;

  SmallVector(void);

  explicit SmallVector(const std::size_t a_initial_size);

  SmallVector(const std::size_t a_initial_size, const T& a_value);

  explicit SmallVector(std::initializer_list<T> init);

  T& operator[](const std::size_t a_index);

  const T& operator[](const std::size_t a_index) const;

  std::size_t size(void) const;

  std::size_t capacity(void) const;

  static constexpr std::size_t staticAllocationCapacity(void);
  void resize(const std::size_t a_count, T val = T());

  void reserve(const std::size_t a_count);

  void assign(const std::size_t a_count, const T& a_value);

  T& front(void);

  const T& front(void) const;

  T& back(void);

  const T& back(void) const;

  bool empty(void) const;
  bool full(void) const;

  void pop_back(void);

  void insert(const_iterator a_pos, const T& a_object);
  void emplace(const_iterator a_pos, T&& a_object);

  void push_back(const T& a_object);
  void emplace_back(T&& a_object);

  void erase(iterator a_pos);
  void erase(iterator a_start, iterator a_end);

  void clear(void);

  // Copy constructor
  SmallVector(const SmallVector& a_rhs);

  // Copy assignment
  SmallVector& operator=(const SmallVector& a_rhs);

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  ~SmallVector(void) = default;

 private:
  typename BaseVector::allocator_type::arena_type memory_pool_m;
  BaseVector vector_m;
};

}  // namespace IRL

#include "irl/data_structures/small_vector.tpp"

#endif  // #ifndef IRL_NO_ABSL

#endif  // SRC_DATA_STRUCTURES_SMALL_VECTOR_H_
