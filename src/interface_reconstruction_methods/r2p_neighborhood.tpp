// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_R2P_NEIGHBORHOOD_TPP_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_R2P_NEIGHBORHOOD_TPP_

namespace IRL {

template <class CellType>
R2PNeighborhood<CellType>::R2PNeighborhood(void)
    : center_cell_surface_area_m(-1.0),
      center_cell_index_m(static_cast<IRL::UnsignedIndex_t>(-1)) {}

template <class CellType>
void R2PNeighborhood<CellType>::addMember(
    const CellType* a_cell,
    const SeparatedMoments<VolumeMoments>* a_volume_moments) {
  assert(a_cell != nullptr);
  assert(a_volume_moments != nullptr);
  collection_m.push_back(CGD(a_cell, a_volume_moments));
}

template <class CellType>
void R2PNeighborhood<CellType>::emptyNeighborhood(void) {
  collection_m.resize(0);
  center_cell_surface_area_m = -1.0;
  center_cell_index_m = static_cast<IRL::UnsignedIndex_t>(-1);
}

template <class CellType>
void R2PNeighborhood<CellType>::setMember(
    const UnsignedIndex_t a_index, const CellType* a_cell,
    const SeparatedMoments<VolumeMoments>* a_volume_moments) {
  assert(a_cell != nullptr);
  assert(a_volume_moments != nullptr);
  this->checkIndex(a_index);
  collection_m[a_index] = CGD(a_cell, a_volume_moments);
}

template <class CellType>
void R2PNeighborhood<CellType>::setCenterOfStencil(
    const UnsignedIndex_t a_index) {
  this->checkIndex(a_index);
  center_cell_index_m = a_index;
}

template <class CellType>
void R2PNeighborhood<CellType>::setSurfaceArea(const double a_surface_area) {
  center_cell_surface_area_m = a_surface_area;
}

template <class CellType>
const CellType& R2PNeighborhood<CellType>::getCenterCell(void) const {
  this->checkCenterStencilSet();
  return this->getCell(center_cell_index_m);
}

template <class CellType>
const SeparatedMoments<VolumeMoments>&
R2PNeighborhood<CellType>::getCenterCellStoredMoments(void) const {
  this->checkCenterStencilSet();
  return this->getStoredMoments(center_cell_index_m);
}

template <class CellType>
const typename R2PNeighborhood<CellType>::CGD::cell_type&
R2PNeighborhood<CellType>::getCell(const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return collection_m.getCell(a_index);
}

template <class CellType>
const SeparatedMoments<VolumeMoments>&
R2PNeighborhood<CellType>::getStoredMoments(
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return collection_m.getStoredMoments(a_index);
}

template <class CellType>
void R2PNeighborhood<CellType>::resize(const UnsignedIndex_t a_size) {
  collection_m.resize(a_size);
}

template <class CellType>
UnsignedIndex_t R2PNeighborhood<CellType>::size(void) const {
  return static_cast<UnsignedIndex_t>(collection_m.size());
}

template <class CellType>
double R2PNeighborhood<CellType>::getSurfaceArea(void) const {
  assert(center_cell_surface_area_m > 0.0);
  return center_cell_surface_area_m;
}

template <class CellType>
typename R2PNeighborhood<CellType>::iterator R2PNeighborhood<CellType>::begin(
    void) noexcept {
  return collection_m.begin();
}
template <class CellType>
typename R2PNeighborhood<CellType>::const_iterator
R2PNeighborhood<CellType>::begin(void) const noexcept {
  return this->cbegin();
}
template <class CellType>
typename R2PNeighborhood<CellType>::const_iterator
R2PNeighborhood<CellType>::end(void) const noexcept {
  return this->cend();
}
template <class CellType>
typename R2PNeighborhood<CellType>::const_iterator
R2PNeighborhood<CellType>::cbegin(void) const noexcept {
  return collection_m.cbegin();
}
template <class CellType>
typename R2PNeighborhood<CellType>::iterator R2PNeighborhood<CellType>::end(
    void) noexcept {
  return collection_m.end();
}
template <class CellType>
typename R2PNeighborhood<CellType>::const_iterator
R2PNeighborhood<CellType>::cend(void) const noexcept {
  return collection_m.cend();
}

template <class CellType>
void R2PNeighborhood<CellType>::checkIndex(UnsignedIndex_t a_index) const {
  assert(a_index < collection_m.size());
}

template <class CellType>
void R2PNeighborhood<CellType>::checkCenterStencilSet(void) const {
  assert(center_cell_index_m != static_cast<UnsignedIndex_t>(-1));
}

}  // namespace IRL

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_R2P_NEIGHBORHOOD_TPP_
