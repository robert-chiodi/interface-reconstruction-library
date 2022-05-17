// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PLVIRA_NEIGHBORHOOD_TPP_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PLVIRA_NEIGHBORHOOD_TPP_

namespace IRL {

template <class CellType>
PLVIRANeighborhood<CellType>::PLVIRANeighborhood(void)
    : center_cell_index_m(static_cast<IRL::UnsignedIndex_t>(-1)) {}

template <class CellType>
void PLVIRANeighborhood<CellType>::addMember(const CellType* a_cell,
                                             const double* a_volume_fraction,
                                             const double a_weight) {
  assert(a_cell != nullptr);
  assert(a_volume_fraction != nullptr);
  collection_m.push_back(CGD(a_cell, a_volume_fraction));
  weights_m.push_back(a_weight);
}

template <class CellType>
void PLVIRANeighborhood<CellType>::emptyNeighborhood(void) {
  collection_m.resize(0);
  center_cell_index_m = static_cast<IRL::UnsignedIndex_t>(-1);
}

template <class CellType>
void PLVIRANeighborhood<CellType>::setCenterOfStencil(
    const UnsignedIndex_t a_index) {
  this->checkIndex(a_index);
  center_cell_index_m = a_index;
}

template <class CellType>
void PLVIRANeighborhood<CellType>::setConstrainedCell(
    CellType* a_cell, const VolumeMoments a_constraints) {
  constrained_cell_m = a_cell;
  constraints_m = a_constraints;
}

template <class CellType>
const CellType& PLVIRANeighborhood<CellType>::getCenterCell(void) const {
  this->checkCenterStencilSet();
  return this->getCell(center_cell_index_m);
}

template <class CellType>
const CellType& PLVIRANeighborhood<CellType>::getConstrainedCell(void) const {
  return *constrained_cell_m;
}

template <class CellType>
const VolumeMoments& PLVIRANeighborhood<CellType>::getConstraints(void) const {
  return constraints_m;
}

template <class CellType>
const double& PLVIRANeighborhood<CellType>::getWeight(
    const UnsignedIndex_t a_index) const {
  return weights_m[a_index];
}

template <class CellType>
UnsignedIndex_t PLVIRANeighborhood<CellType>::getCenterOfStencilIndex(
    void) const {
  this->checkCenterStencilSet();
  return center_cell_index_m;
}

template <class CellType>
const double& PLVIRANeighborhood<CellType>::getCenterCellStoredMoments(
    void) const {
  this->checkCenterStencilSet();
  return this->getStoredMoments(center_cell_index_m);
}

template <class CellType>
const typename PLVIRANeighborhood<CellType>::CGD::cell_type&
PLVIRANeighborhood<CellType>::getCell(const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return collection_m.getCell(a_index);
}

template <class CellType>
const double& PLVIRANeighborhood<CellType>::getStoredMoments(
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return collection_m.getStoredMoments(a_index);
}

template <class CellType>
void PLVIRANeighborhood<CellType>::resize(const UnsignedIndex_t a_size) {
  collection_m.resize(a_size);
}

template <class CellType>
UnsignedIndex_t PLVIRANeighborhood<CellType>::size(void) const {
  return static_cast<UnsignedIndex_t>(collection_m.size());
}

template <class CellType>
typename PLVIRANeighborhood<CellType>::iterator
PLVIRANeighborhood<CellType>::begin(void) noexcept {
  return collection_m.begin();
}

template <class CellType>
typename PLVIRANeighborhood<CellType>::const_iterator
PLVIRANeighborhood<CellType>::begin(void) const noexcept {
  return this->cbegin();
}

template <class CellType>
typename PLVIRANeighborhood<CellType>::const_iterator
PLVIRANeighborhood<CellType>::end(void) const noexcept {
  return this->cend();
}

template <class CellType>
typename PLVIRANeighborhood<CellType>::const_iterator
PLVIRANeighborhood<CellType>::cbegin(void) const noexcept {
  return collection_m.cbegin();
}

template <class CellType>
typename PLVIRANeighborhood<CellType>::iterator
PLVIRANeighborhood<CellType>::end(void) noexcept {
  return collection_m.end();
}

template <class CellType>
typename PLVIRANeighborhood<CellType>::const_iterator
PLVIRANeighborhood<CellType>::cend(void) const noexcept {
  return collection_m.cend();
}

template <class CellType>
void PLVIRANeighborhood<CellType>::checkIndex(UnsignedIndex_t a_index) const {
  assert(a_index < collection_m.size());
}

template <class CellType>
void PLVIRANeighborhood<CellType>::checkCenterStencilSet(void) const {
  assert(center_cell_index_m != static_cast<UnsignedIndex_t>(-1));
}

}  // namespace IRL

#endif  // IRL_INTERFACE_RECONSTRUCTION_METHODS_PLVIRA_NEIGHBORHOOD_TPP_
