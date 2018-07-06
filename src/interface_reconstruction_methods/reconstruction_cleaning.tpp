// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_CLEANING_TPP_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_CLEANING_TPP_

namespace IRL {

template <class CellType, class ReconstructionTyoe>
inline void cleanReconstruction(const CellType& a_cell,
                                const double a_volume_fraction,
                                ReconstructionTyoe* a_reconstruction) {
  cleanReconstructionOutOfCell(a_cell, a_volume_fraction, a_reconstruction);
  cleanReconstructionSameNormal(a_cell, a_volume_fraction, a_reconstruction);
}

template <class CellType, class ReconstructionType>
inline void cleanReconstructionOutOfCell(const CellType& a_cell,
                                         const double a_volume_fraction,
                                         ReconstructionType* a_reconstruction) {
  for (UnsignedIndex_t n = a_reconstruction->getNumberOfPlanes() - 1;
       n != static_cast<UnsignedIndex_t>(-1); --n) {
    if (!isPlaneIntersectingCell((*a_reconstruction)[n], a_cell)) {
      a_reconstruction->removePlane(n);
    }
  }
  // If all planes removed, set one plane that recreates correct phase.
  if (a_reconstruction->getNumberOfPlanes() == 0) {
    setToPurePhaseReconstruction(a_volume_fraction, a_reconstruction);
  }
}

template <class CellType, class ReconstructionType>
inline void cleanReconstructionSameNormal(
    const CellType& a_cell, const double a_volume_fraction,
    ReconstructionType* a_reconstruction) {
  auto starting_number_of_planes = a_reconstruction->getNumberOfPlanes();
  // Check if two normals but they're the same
  auto it = std::unique(a_reconstruction->begin(), a_reconstruction->end());
  a_reconstruction->setNumberOfPlanes(
      std::max(static_cast<UnsignedIndex_t>(1),
               static_cast<UnsignedIndex_t>(
                   std::distance(a_reconstruction->begin(), it))));
  if (a_reconstruction->getNumberOfPlanes() < starting_number_of_planes) {
    setDistanceToMatchVolumeFraction(a_cell, a_volume_fraction,
                                     a_reconstruction);
  }
}

template <class ReconstructionType>
void setToPurePhaseReconstruction(const double a_internal_volume_fraction,
                                  ReconstructionType* a_reconstruction) {
  *a_reconstruction = ReconstructionType::fromOnePlane(
      Plane(Normal(0.0, 0.0, 0.0),
            std::copysign(global_constants::ARBITRARILY_LARGE_DISTANCE,
                          a_internal_volume_fraction - 0.5)));
}

}  // namespace IRL

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_CLEANING_TPP_
