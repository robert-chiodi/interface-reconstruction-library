// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_STACK_VECTOR_TPP_
#define SRC_DATA_STRUCTURES_STACK_VECTOR_TPP_

#include <algorithm>

namespace IRL {

template <class ObjectType, UnsignedIndex_t kMaxSize>
StackVector<ObjectType, kMaxSize>::StackVector(void) : size_m(0) {}

template <class ObjectType, UnsignedIndex_t kMaxSize>
template <UnsignedIndex_t kOtherMaxSize>
StackVector<ObjectType, kMaxSize>::StackVector(
    const StackVector<ObjectType, kOtherMaxSize>& other)
    : size_m(other.size()) {
  assert(other.size() <= this->maxSize());
  std::copy(other.begin(), other.end(), this->begin());
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
template <UnsignedIndex_t kOtherMaxSize>
StackVector<ObjectType, kMaxSize>& StackVector<ObjectType, kMaxSize>::operator=(
    const StackVector<ObjectType, kOtherMaxSize>& other) {
  assert(other.size() <= this->maxSize());
  this->resize(other.size());
  std::copy(other.begin(), other.end(), this->begin());
  return (*this);
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
ObjectType& StackVector<ObjectType, kMaxSize>::operator[](
    const UnsignedIndex_t a_index) {
  assert(a_index < this->size());
  return object_container_m[a_index];
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
const ObjectType& StackVector<ObjectType, kMaxSize>::operator[](
    const UnsignedIndex_t a_index) const {
  assert(a_index < this->size());
  return object_container_m[a_index];
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
UnsignedIndex_t StackVector<ObjectType, kMaxSize>::size(void) const {
  return size_m;
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
UnsignedIndex_t StackVector<ObjectType, kMaxSize>::maxSize(void) const {
  return static_cast<UnsignedIndex_t>(object_container_m.size());
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
void StackVector<ObjectType, kMaxSize>::resize(const UnsignedIndex_t a_count) {
  assert(a_count <= this->maxSize());
  size_m = a_count;
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
void StackVector<ObjectType, kMaxSize>::assign(const UnsignedIndex_t a_count,
                                               const ObjectType& a_value) {
  this->resize(a_count);
  std::fill(this->begin(), this->end(), a_value);
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
ObjectType& StackVector<ObjectType, kMaxSize>::front(void) {
  assert(!empty());
  return object_container_m[0];
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
const ObjectType& StackVector<ObjectType, kMaxSize>::front(void) const {
  assert(!empty());
  return object_container_m[0];
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
ObjectType& StackVector<ObjectType, kMaxSize>::back(void) {
  assert(!empty());
  return object_container_m[this->size() - 1];
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
const ObjectType& StackVector<ObjectType, kMaxSize>::back(void) const {
  assert(!empty());
  return object_container_m[this->size() - 1];
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
bool StackVector<ObjectType, kMaxSize>::empty(void) const {
  return size_m == 0;
}
template <class ObjectType, UnsignedIndex_t kMaxSize>
bool StackVector<ObjectType, kMaxSize>::full(void) const {
  return size_m == object_container_m.size();
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
void StackVector<ObjectType, kMaxSize>::pop_back(void) {
  assert(!empty());
  this->decrementSize();
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
void StackVector<ObjectType, kMaxSize>::push_back(const ObjectType& a_object) {
  assert(!full());
  object_container_m[this->size()] = a_object;
  this->incrementSize();
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
void StackVector<ObjectType, kMaxSize>::erase(iterator a_pos) {
  assert(this->end() - a_pos > 0);
  while (a_pos != (this->end() - 1)) {
    *a_pos = *(a_pos + 1);
    ++a_pos;
  }
  this->decrementSize();
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
typename StackVector<ObjectType, kMaxSize>::iterator
StackVector<ObjectType, kMaxSize>::begin(void) noexcept {
  return object_container_m.begin();
}
template <class ObjectType, UnsignedIndex_t kMaxSize>
typename StackVector<ObjectType, kMaxSize>::const_iterator
StackVector<ObjectType, kMaxSize>::begin(void) const noexcept {
  return this->cbegin();
}
template <class ObjectType, UnsignedIndex_t kMaxSize>
typename StackVector<ObjectType, kMaxSize>::const_iterator
StackVector<ObjectType, kMaxSize>::cbegin(void) const noexcept {
  return object_container_m.cbegin();
}
template <class ObjectType, UnsignedIndex_t kMaxSize>
typename StackVector<ObjectType, kMaxSize>::iterator
StackVector<ObjectType, kMaxSize>::end(void) noexcept {
  return this->begin() + this->size();
}
template <class ObjectType, UnsignedIndex_t kMaxSize>
typename StackVector<ObjectType, kMaxSize>::const_iterator
StackVector<ObjectType, kMaxSize>::end(void) const noexcept {
  return this->cend();
}
template <class ObjectType, UnsignedIndex_t kMaxSize>
typename StackVector<ObjectType, kMaxSize>::const_iterator
StackVector<ObjectType, kMaxSize>::cend(void) const noexcept {
  return object_container_m.cbegin() + this->size();
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
void StackVector<ObjectType, kMaxSize>::incrementSize(void) {
  assert(!full());
  ++size_m;
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
void StackVector<ObjectType, kMaxSize>::decrementSize(void) {
  assert(!empty());
  --size_m;
}

template <class ObjectType, UnsignedIndex_t kMaxSize>
inline std::ostream& operator<<(
    std::ostream& out,
    const StackVector<ObjectType, kMaxSize>& a_stack_vector) {
  out << "There are " << a_stack_vector.size()
      << " elements in the StackVector. " << std::endl;
  out << std::setprecision(15);
  for (const auto& element : a_stack_vector) {
    out << element << " ";
  }
  out << std::endl;
  return out;
}

}  // namespace IRL

#endif  // SRC_DATA_STRUCTURES_STACK_VECTOR_TPP_
