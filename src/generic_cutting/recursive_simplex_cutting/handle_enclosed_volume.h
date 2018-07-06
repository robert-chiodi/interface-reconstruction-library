// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_HANDLE_ENCLOSED_VOLUME_H_
#define SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_HANDLE_ENCLOSED_VOLUME_H_

#include "src/generic_cutting/generic_cutting_definitions.h"

namespace IRL {

struct AccumulateIntoScalar {
  template <class ReturnType, class ReconstructionType, class SimplexType>
  inline static void handleVolumeCompletelyInternalToPlanarReconstruction(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

struct PassToNestedType {
  template <class ReturnType, class ReconstructionType, class SimplexType>
  inline static void handleVolumeCompletelyInternalToPlanarReconstruction(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

struct AccumulateIntoCollection {
  template <class ReturnType, class ReconstructionType, class SimplexType>
  inline static void handleVolumeCompletelyInternalToPlanarReconstruction(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

}  // namespace IRL

#include "src/generic_cutting/recursive_simplex_cutting/handle_enclosed_volume.tpp"

#endif  // SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_HANDLE_ENCLOSED_VOLUME_H_
