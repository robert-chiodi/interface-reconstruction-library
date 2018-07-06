// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_LVIRA_NEIGHBORHOOD_TPP_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_LVIRA_NEIGHBORHOOD_TPP_

namespace IRL {

template <class CellType>
LVIRANeighborhood<CellType>::LVIRANeighborhood(void)
    : center_cell_index_m(static_cast<IRL::UnsignedIndex_t>(-1)) {}

template <class CellType>
void LVIRANeighborhood<CellType>::addMember(const CellType* a_cell,
                                            const double* a_volume_fraction) {
  assert(a_cell != nullptr);
  assert(a_volume_fraction != nullptr);
  collection_m.push_back(CGD(a_cell, a_volume_fraction));
}

template <class CellType>
void LVIRANeighborhood<CellType>::emptyNeighborhood(void) {
  collection_m.resize(0);
  center_cell_index_m = static_cast<IRL::UnsignedIndex_t>(-1);
}

template <class CellType>
void LVIRANeighborhood<CellType>::setMember(const UnsignedIndex_t a_index,
                                            const CellType* a_cell,
                                            const double* a_volume_fraction) {
  assert(a_cell != nullptr);
  assert(a_volume_fraction != nullptr);
  this->checkIndex(a_index);
  collection_m[a_index] = CGD(a_cell, a_volume_fraction);
}

template <class CellType>
void LVIRANeighborhood<CellType>::setCenterOfStencil(
    const UnsignedIndex_t a_index) {
  this->checkIndex(a_index);
  center_cell_index_m = a_index;
}

template <class CellType>
const CellType& LVIRANeighborhood<CellType>::getCenterCell(void) const {
  this->checkCenterStencilSet();
  return this->getCell(center_cell_index_m);
}

template <class CellType>
UnsignedIndex_t LVIRANeighborhood<CellType>::getCenterOfStencilIndex(
    void) const {
  this->checkCenterStencilSet();
  return center_cell_index_m;
}

template <class CellType>
const double& LVIRANeighborhood<CellType>::getCenterCellStoredMoments(
    void) const {
  this->checkCenterStencilSet();
  return this->getStoredMoments(center_cell_index_m);
}

template <class CellType>
const typename LVIRANeighborhood<CellType>::CGD::cell_type&
LVIRANeighborhood<CellType>::getCell(const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return collection_m.getCell(a_index);
}

template <class CellType>
const double& LVIRANeighborhood<CellType>::getStoredMoments(
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return collection_m.getStoredMoments(a_index);
}

template <class CellType>
void LVIRANeighborhood<CellType>::resize(const UnsignedIndex_t a_size) {
  collection_m.resize(a_size);
}

template <class CellType>
UnsignedIndex_t LVIRANeighborhood<CellType>::size(void) const {
  return static_cast<UnsignedIndex_t>(collection_m.size());
}

template <class CellType>
typename LVIRANeighborhood<CellType>::iterator
LVIRANeighborhood<CellType>::begin(void) noexcept {
  return collection_m.begin();
}

template <class CellType>
typename LVIRANeighborhood<CellType>::const_iterator
LVIRANeighborhood<CellType>::begin(void) const noexcept {
  return this->cbegin();
}

template <class CellType>
typename LVIRANeighborhood<CellType>::const_iterator
LVIRANeighborhood<CellType>::end(void) const noexcept {
  return this->cend();
}

template <class CellType>
typename LVIRANeighborhood<CellType>::const_iterator
LVIRANeighborhood<CellType>::cbegin(void) const noexcept {
  return collection_m.cbegin();
}

template <class CellType>
typename LVIRANeighborhood<CellType>::iterator LVIRANeighborhood<CellType>::end(
    void) noexcept {
  return collection_m.end();
}

template <class CellType>
typename LVIRANeighborhood<CellType>::const_iterator
LVIRANeighborhood<CellType>::cend(void) const noexcept {
  return collection_m.cend();
}

template <class CellType>
void LVIRANeighborhood<CellType>::checkIndex(UnsignedIndex_t a_index) const {
  assert(a_index < collection_m.size());
}

template <class CellType>
void LVIRANeighborhood<CellType>::checkCenterStencilSet(void) const {
  assert(center_cell_index_m != static_cast<UnsignedIndex_t>(-1));
}

}  // namespace IRL

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_LVIRA_NEIGHBORHOOD_TPP_
