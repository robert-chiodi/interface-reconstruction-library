// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_CELL_COLLECTION_TPP_
#define SRC_MOMENTS_CELL_COLLECTION_TPP_

namespace IRL {

template <class CellGroupType>
const typename CellGroupType::cell_type& CellCollection<CellGroupType>::getCell(
    const UnsignedIndex_t a_index) const {
  return (*this)[a_index].getCell();
}

template <class CellGroupType>
const typename CellGroupType::contained_type&
CellCollection<CellGroupType>::getStoredMoments(
    const UnsignedIndex_t a_index) const {
  return (*this)[a_index].getStoredMoments();
}

template <class CellGroupType>
UnsignedIndex_t CellCollection<CellGroupType>::size(void) const {
  return collection_m.size();
}

template <class CellGroupType>
void CellCollection<CellGroupType>::resize(const UnsignedIndex_t a_size) {
  collection_m.resize(a_size);
}

template <class CellGroupType>
void CellCollection<CellGroupType>::push_back(const CellGroupType& a_object) {
  collection_m.push_back(a_object);
}

template <class CellGroupType>
void CellCollection<CellGroupType>::clear(void) {
  collection_m.clear();
}

template <class CellGroupType>
CellGroupType& CellCollection<CellGroupType>::operator[](
    const UnsignedIndex_t a_index) {
  this->checkIndex(a_index);
  return collection_m[a_index];
}

template <class CellGroupType>
const CellGroupType& CellCollection<CellGroupType>::operator[](
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return collection_m[a_index];
}

template <class CellGroupType>
typename CellCollection<CellGroupType>::iterator
CellCollection<CellGroupType>::begin(void) noexcept {
  return collection_m.begin();
}
template <class CellGroupType>
typename CellCollection<CellGroupType>::const_iterator
CellCollection<CellGroupType>::begin(void) const noexcept {
  return this->cbegin();
}
template <class CellGroupType>
typename CellCollection<CellGroupType>::const_iterator
CellCollection<CellGroupType>::cbegin(void) const noexcept {
  return collection_m.cbegin();
}
template <class CellGroupType>
typename CellCollection<CellGroupType>::iterator
CellCollection<CellGroupType>::end(void) noexcept {
  return collection_m.end();
}
template <class CellGroupType>
typename CellCollection<CellGroupType>::const_iterator
CellCollection<CellGroupType>::end(void) const noexcept {
  return this->cend();
}
template <class CellGroupType>
typename CellCollection<CellGroupType>::const_iterator
CellCollection<CellGroupType>::cend(void) const noexcept {
  return collection_m.cend();
}

template <class CellGroupType>
void CellCollection<CellGroupType>::checkIndex(
    const UnsignedIndex_t a_index) const {
  assert(a_index < this->size());
}

}  // namespace IRL

#endif  // SRC_MOMENTS_CELL_COLLECTION_TPP_
