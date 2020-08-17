// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_COLLECTION_TPP_
#define SRC_DATA_STRUCTURES_COLLECTION_TPP_

namespace IRL {

template <class ObjectType>
UnsignedIndex_t Collection<ObjectType>::size(void) const {
  return static_cast<UnsignedIndex_t>(collection_m.size());
}

template <class ObjectType>
void Collection<ObjectType>::resize(const UnsignedIndex_t a_size) {
  collection_m.resize(a_size);
}

template <class ObjectType>
void Collection<ObjectType>::push_back(const ObjectType& a_object) {
  collection_m.push_back(a_object);
}

template <class ObjectType>
void Collection<ObjectType>::emplace_back(const ObjectType& a_object) {
  collection_m.emplace_back(a_object);
}

template <class ObjectType>
void Collection<ObjectType>::clear(void) {
  collection_m.clear();
}

template <class ObjectType>
void Collection<ObjectType>::erase(const UnsignedIndex_t a_index) {
  collection_m.erase(this->begin() + a_index);
}

template <class ObjectType>
ObjectType& Collection<ObjectType>::operator[](const UnsignedIndex_t a_index) {
  this->checkIndex(a_index);
  return collection_m[a_index];
}

template <class ObjectType>
const ObjectType& Collection<ObjectType>::operator[](
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return collection_m[a_index];
}

template <class ObjectType>
typename Collection<ObjectType>::iterator Collection<ObjectType>::begin(
    void) noexcept {
  return collection_m.begin();
}
template <class ObjectType>
typename Collection<ObjectType>::const_iterator Collection<ObjectType>::begin(
    void) const noexcept {
  return this->cbegin();
}
template <class ObjectType>
typename Collection<ObjectType>::const_iterator Collection<ObjectType>::cbegin(
    void) const noexcept {
  return collection_m.cbegin();
}
template <class ObjectType>
typename Collection<ObjectType>::iterator Collection<ObjectType>::end(
    void) noexcept {
  return collection_m.end();
}
template <class ObjectType>
typename Collection<ObjectType>::const_iterator Collection<ObjectType>::end(
    void) const noexcept {
  return this->cend();
}
template <class ObjectType>
typename Collection<ObjectType>::const_iterator Collection<ObjectType>::cend(
    void) const noexcept {
  return collection_m.cend();
}

template <class ObjectType>
void Collection<ObjectType>::checkIndex(const UnsignedIndex_t a_index) const {
  assert(a_index < this->size());
}

}  // namespace IRL

#endif  // SRC_DATA_STRUCTURES_COLLECTION_TPP_
