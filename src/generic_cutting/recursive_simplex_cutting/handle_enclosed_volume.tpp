// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_HANDLE_ENCLOSED_VOLUME_TPP_
#define SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_HANDLE_ENCLOSED_VOLUME_TPP_

namespace IRL {

template <class ReturnType, class ReconstructionType, class SimplexType>
void AccumulateIntoScalar::handleVolumeCompletelyInternalToPlanarReconstruction(
    const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
    ReturnType* a_moments_to_return) {
  *a_moments_to_return +=
      IRL::getVolumeMoments<ReturnType, RecursiveSimplexCutting>(
          a_simplex, a_reconstruction.getNextReconstruction());
}

template <class ReturnType, class ReconstructionType, class SimplexType>
void PassToNestedType::handleVolumeCompletelyInternalToPlanarReconstruction(
    const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
    ReturnType* a_moments_to_return) {
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  *a_moments_to_return +=
      IRL::getVolumeMoments<WantedVolumeMomentsType, RecursiveSimplexCutting>(
          a_simplex, a_reconstruction.getNextReconstruction());
}

template <class ReturnType, class ReconstructionType, class SimplexType>
void AccumulateIntoCollection::
    handleVolumeCompletelyInternalToPlanarReconstruction(
        const SimplexType& a_simplex,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  assert(a_reconstruction.isIdSet());
  (*a_moments_to_return)[a_reconstruction.getId()] +=
      IRL::getVolumeMoments<WantedVolumeMomentsType, RecursiveSimplexCutting>(
          a_simplex, a_reconstruction.getNextReconstruction());
}

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_HANDLE_ENCLOSED_VOLUME_TPP_
