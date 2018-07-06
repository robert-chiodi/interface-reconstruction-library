// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/interface_reconstruction_methods/elvira_neighborhood.h"

namespace IRL {

/// \brief Construct a CellGroupedMoments and add it
/// to the collection for index i,j,k.
void ELVIRANeighborhood::setMember(
    const RectangularCuboid* a_rectangular_cuboid,
    const double* a_liquid_volume_fraction, const int i, const int j,
    const int k) {
  assert(a_rectangular_cuboid != nullptr);
  assert(a_liquid_volume_fraction != nullptr);
  collection_m[this->calculateLinearIndex(i, j, k)] =
      CGD(a_rectangular_cuboid, a_liquid_volume_fraction);
}

/// \brief Return the cell stored at the index i,j,k
const RectangularCuboid& ELVIRANeighborhood::getCell(const int i, const int j,
                                                     const int k) const {
  return collection_m.getCell(this->calculateLinearIndex(i, j, k));
}

/// \brief Return moments stored at the index i,j,k
double ELVIRANeighborhood::getStoredMoments(const int i, const int j,
                                            const int k) const {
  return collection_m.getStoredMoments(this->calculateLinearIndex(i, j, k));
}

/// \brief Set size of the neighborhood.
void ELVIRANeighborhood::resize(const UnsignedIndex_t a_size) {
  collection_m.resize(a_size);
}

/// \brief Calculate linear index from i,j,k.
UnsignedIndex_t ELVIRANeighborhood::calculateLinearIndex(const int i,
                                                         const int j,
                                                         const int k) const {
  return static_cast<UnsignedIndex_t>((i + 1) + (j + 1) * 3 + (k + 1) * 9);
}

}  // namespace IRL
