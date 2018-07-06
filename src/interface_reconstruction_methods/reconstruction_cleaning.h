// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_CLEANING_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_CLEANING_H_

#include <algorithm>

#include "src/interface_reconstruction_methods/volume_fraction_matching.h"
#include "src/parameters/constants.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class CellType, class ReconstructionTyoe>
inline void cleanReconstruction(const CellType& a_cell,
                                const double a_volume_fraction,
                                ReconstructionTyoe* a_reconstruction);

template <class CellType, class ReconstructionType>
inline void cleanReconstructionOutOfCell(const CellType& a_cell,
                                         const double a_volume_fraction,
                                         ReconstructionType* a_reconstruction);

template <class CellType, class ReconstructionType>
inline void cleanReconstructionSameNormal(const CellType& a_cell,
                                          const double a_volume_fraction,
                                          ReconstructionType* a_reconstruction);

template <class ReconstructionType>
inline void setToPurePhaseReconstruction(
    const double a_internal_volume_fraction,
    ReconstructionType* a_reconstruction);

}  // namespace IRL

#include "src/interface_reconstruction_methods/reconstruction_cleaning.tpp"

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_RECONSTRUCTION_CLEANING_H_
