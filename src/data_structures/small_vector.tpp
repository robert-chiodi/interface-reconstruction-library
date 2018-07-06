// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_SMALL_VECTOR_TPP_
#define SRC_DATA_STRUCTURES_SMALL_VECTOR_TPP_

namespace IRL {

template <class T, UnsignedIndex_t BuffSizeInElements>
SmallVector<T, BuffSizeInElements>::SmallVector(void)
    : memory_pool_m(), vector_m(memory_pool_m) {
  this->reserve(BuffSizeInElements);
}

template <class T, UnsignedIndex_t BuffSizeInElements>
SmallVector<T, BuffSizeInElements>::SmallVector(
    const UnsignedIndex_t a_initial_size)
    : memory_pool_m(), vector_m(a_initial_size, T(), memory_pool_m) {
  this->reserve(BuffSizeInElements);
}

template <class T, UnsignedIndex_t BuffSizeInElements>
SmallVector<T, BuffSizeInElements>::SmallVector(
    const UnsignedIndex_t a_initial_size, const T& a_value)
    : memory_pool_m(), vector_m(a_initial_size, a_value, memory_pool_m) {
  this->reserve(BuffSizeInElements);
}

template <class T, UnsignedIndex_t BuffSizeInElements>
SmallVector<T, BuffSizeInElements>::SmallVector(std::initializer_list<T> init)
    : memory_pool_m(), vector_m(init, memory_pool_m) {
  this->reserve(BuffSizeInElements);
}
template <class T, UnsignedIndex_t BuffSizeInElements>
T& SmallVector<T, BuffSizeInElements>::operator[](
    const UnsignedIndex_t a_index) {
  return vector_m[a_index];
}

template <class T, UnsignedIndex_t BuffSizeInElements>
const T& SmallVector<T, BuffSizeInElements>::operator[](
    const UnsignedIndex_t a_index) const {
  return vector_m[a_index];
}

template <class T, UnsignedIndex_t BuffSizeInElements>
UnsignedIndex_t SmallVector<T, BuffSizeInElements>::size(void) const {
  return static_cast<UnsignedIndex_t>(vector_m.size());
}

template <class T, UnsignedIndex_t BuffSizeInElements>
UnsignedIndex_t SmallVector<T, BuffSizeInElements>::capacity(void) const {
  return static_cast<UnsignedIndex_t>(vector_m.capacity());
}

template <class T, UnsignedIndex_t BuffSizeInElements>
constexpr UnsignedIndex_t
SmallVector<T, BuffSizeInElements>::staticAllocationCapacity(void) {
  return BuffSizeInElements;
}

template <class T, UnsignedIndex_t BuffSizeInElements>
void SmallVector<T, BuffSizeInElements>::resize(const UnsignedIndex_t a_count,
                                                T val) {
  vector_m.resize(a_count, val);
}

template <class T, UnsignedIndex_t BuffSizeInElements>
void SmallVector<T, BuffSizeInElements>::reserve(
    const UnsignedIndex_t a_count) {
  vector_m.reserve(a_count);
}

template <class T, UnsignedIndex_t BuffSizeInElements>
void SmallVector<T, BuffSizeInElements>::assign(const UnsignedIndex_t a_count,
                                                const T& a_value) {
  vector_m.assign(a_count, a_value);
}

template <class T, UnsignedIndex_t BuffSizeInElements>
T& SmallVector<T, BuffSizeInElements>::front(void) {
  return vector_m;
}

template <class T, UnsignedIndex_t BuffSizeInElements>
const T& SmallVector<T, BuffSizeInElements>::front(void) const {
  return vector_m;
}

template <class T, UnsignedIndex_t BuffSizeInElements>
T& SmallVector<T, BuffSizeInElements>::back(void) {
  return vector_m.front();
}

template <class T, UnsignedIndex_t BuffSizeInElements>
const T& SmallVector<T, BuffSizeInElements>::back(void) const {
  return vector_m.back();
}

template <class T, UnsignedIndex_t BuffSizeInElements>
bool SmallVector<T, BuffSizeInElements>::empty(void) const {
  return vector_m.empty();
}
template <class T, UnsignedIndex_t BuffSizeInElements>
bool SmallVector<T, BuffSizeInElements>::full(void) const {
  return vector_m.full();
}

template <class T, UnsignedIndex_t BuffSizeInElements>
void SmallVector<T, BuffSizeInElements>::pop_back(void) {
  vector_m.pop_back();
}

template <class T, UnsignedIndex_t BuffSizeInElements>
void SmallVector<T, BuffSizeInElements>::insert(const_iterator a_pos,
                                                const T& a_object) {
  vector_m.insert(a_pos, a_object);
}
template <class T, UnsignedIndex_t BuffSizeInElements>
void SmallVector<T, BuffSizeInElements>::emplace(const_iterator a_pos,
                                                 T&& a_object) {
  vector_m.emaplce(a_pos, a_object);
}

template <class T, UnsignedIndex_t BuffSizeInElements>
void SmallVector<T, BuffSizeInElements>::push_back(const T& a_object) {
  vector_m.push_back(a_object);
}
template <class T, UnsignedIndex_t BuffSizeInElements>
void SmallVector<T, BuffSizeInElements>::emplace_back(T&& a_object) {
  vector_m.emplace_back(std::move(a_object));
}

template <class T, UnsignedIndex_t BuffSizeInElements>
void SmallVector<T, BuffSizeInElements>::erase(iterator a_pos) {
  vector_m.erase(a_pos);
}
template <class T, UnsignedIndex_t BuffSizeInElements>
void SmallVector<T, BuffSizeInElements>::erase(iterator a_start,
                                               iterator a_end) {
  vector_m.erase(a_start, a_end);
}

template <class T, UnsignedIndex_t BuffSizeInElements>
void SmallVector<T, BuffSizeInElements>::clear(void) {
  vector_m.clear();
}

template <class T, UnsignedIndex_t BuffSizeInElements>
SmallVector<T, BuffSizeInElements>::SmallVector(const SmallVector& a_rhs)
    : memory_pool_m(), vector_m{memory_pool_m} {
  vector_m = a_rhs.vector_m;
}

template <class T, UnsignedIndex_t BuffSizeInElements>
SmallVector<T, BuffSizeInElements>& SmallVector<T, BuffSizeInElements>::
operator=(const SmallVector& a_rhs) {
  if (this != &a_rhs) {
    vector_m = a_rhs.vector_m;
  }
  return (*this);
}

template <class T, UnsignedIndex_t BuffSizeInElements>
typename SmallVector<T, BuffSizeInElements>::iterator
SmallVector<T, BuffSizeInElements>::begin(void) noexcept {
  return vector_m.begin();
}
template <class T, UnsignedIndex_t BuffSizeInElements>
typename SmallVector<T, BuffSizeInElements>::const_iterator
SmallVector<T, BuffSizeInElements>::begin(void) const noexcept {
  return this->cbegin();
}
template <class T, UnsignedIndex_t BuffSizeInElements>
typename SmallVector<T, BuffSizeInElements>::const_iterator
SmallVector<T, BuffSizeInElements>::cbegin(void) const noexcept {
  return vector_m.cbegin();
}
template <class T, UnsignedIndex_t BuffSizeInElements>
typename SmallVector<T, BuffSizeInElements>::iterator
SmallVector<T, BuffSizeInElements>::end(void) noexcept {
  return vector_m.end();
}
template <class T, UnsignedIndex_t BuffSizeInElements>
typename SmallVector<T, BuffSizeInElements>::const_iterator
SmallVector<T, BuffSizeInElements>::end(void) const noexcept {
  return this->cend();
}
template <class T, UnsignedIndex_t BuffSizeInElements>
typename SmallVector<T, BuffSizeInElements>::const_iterator
SmallVector<T, BuffSizeInElements>::cend(void) const noexcept {
  return vector_m.cend();
}

}  // namespace IRL

#endif  // SRC_DATA_STRUCTURES_SMALL_VECTOR_TPP_
