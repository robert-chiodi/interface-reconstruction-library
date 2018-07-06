// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_CELL_GROUPED_MOMENTS_TPP_
#define SRC_MOMENTS_CELL_GROUPED_MOMENTS_TPP_

namespace IRL {

template <class CellType, class ContainedMomentsType>
CellGroupedMoments<CellType, ContainedMomentsType>::CellGroupedMoments(void)
    : ptr_to_cell_m(nullptr), ptr_to_contained_m(nullptr) {}

template <class CellType, class ContainedMomentsType>
CellGroupedMoments<CellType, ContainedMomentsType>::CellGroupedMoments(
    const CellType* a_cell, const ContainedMomentsType* a_contained_object)
    : ptr_to_cell_m(a_cell), ptr_to_contained_m(a_contained_object) {
  assert(a_cell != nullptr);
  assert(a_contained_object != nullptr);
}

template <class CellType, class ContainedMomentsType>
const CellType& CellGroupedMoments<CellType, ContainedMomentsType>::getCell(
    void) const {
  assert(ptr_to_cell_m != nullptr);
  return *ptr_to_cell_m;
}

template <class CellType, class ContainedMomentsType>
const ContainedMomentsType&
CellGroupedMoments<CellType, ContainedMomentsType>::getStoredMoments(
    void) const {
  assert(ptr_to_contained_m != nullptr);
  return *ptr_to_contained_m;
}

template <class CellType, class ContainedMomentsType>
ContainedMomentsType CellGroupedMoments<CellType, ContainedMomentsType>::
    calculateNormalizedVolumeMoments(const PlanarSeparator& a_separator) const {
  return getNormalizedVolumeMoments<ContainedMomentsType,
                                    ReconstructionDefaultCuttingMethod>(
      this->getCell(), a_separator);
}

}  // namespace IRL

#endif  // SRC_MOMENTS_CELL_GROUPED_MOMENTS_TPP_
