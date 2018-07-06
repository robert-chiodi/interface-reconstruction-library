// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CUT_SIMPLEX_DRIVERS_H_
#define SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CUT_SIMPLEX_DRIVERS_H_

#include "src/generic_cutting/recursive_simplex_cutting/continue_dividing_volume.h"
#include "src/generic_cutting/recursive_simplex_cutting/handle_enclosed_volume.h"
#include "src/parameters/defined_types.h"

namespace IRL {

// Main function that drives dividing up a simplex and calculating moments
// for reconstructions. Allows changing handling of the volune internal to the
// reconstruction through choosing a
// HowToHandleVolumeCompletelyInternalToPlanarReconstruction structure from
// src/generic_cutting/handle_enclosed_volume.h and how to carry out division of
// a simplex according to HowToContinueDividingVolumeByPlanarReconstruction
// chosen from src/generic_cutting/continue_dividing_volume.h.
template <class HowToHandleVolumeCompletelyInternalToPlanarReconstruction,
          class HowToContinueDividingVolumeByPlanarReconstruction>
class DivideSimplexByPlanarReconstruction
    : public HowToHandleVolumeCompletelyInternalToPlanarReconstruction,
      public HowToContinueDividingVolumeByPlanarReconstruction {
 public:
  template <class SimplexType, class ReconstructionType, class ReturnType>
  static void execute(const SimplexType& a_simplex,
                      const ReconstructionType& a_reconstruction,
                      const UnsignedIndex_t a_cutting_plane_index,
                      ReturnType* a_moments_to_return) {
    if (noPlanesLeftToDivideBy(a_cutting_plane_index, a_reconstruction)) {
      handleVolumeCompletelyInternalToPlanarReconstruction(
          a_simplex, a_reconstruction, a_moments_to_return);
    } else {
      continueDividingVolumeByPlanarReconstruction(a_simplex, a_reconstruction,
                                                   a_cutting_plane_index,
                                                   a_moments_to_return);
    }
  }

 private:
  template <class ReconstructionType>
  static bool noPlanesLeftToDivideBy(
      const UnsignedIndex_t a_cutting_plane_index,
      const ReconstructionType& a_reconstruction) {
    const auto& cutting_reconstruction =
        a_reconstruction.getCurrentReconstruction();

    return a_cutting_plane_index == cutting_reconstruction.getNumberOfPlanes();
  }

  template <class SimplexType, class ReconstructionType, class ReturnType>
  static void handleVolumeCompletelyInternalToPlanarReconstruction(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return) {
    HowToHandleVolumeCompletelyInternalToPlanarReconstruction::
        handleVolumeCompletelyInternalToPlanarReconstruction(
            a_simplex, a_reconstruction, a_moments_to_return);
  }

  template <class SimplexType, class ReconstructionType, class ReturnType>
  static void continueDividingVolumeByPlanarReconstruction(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      ReturnType* a_moments_to_return) {
    HowToContinueDividingVolumeByPlanarReconstruction::
        continueDividingVolumeByPlanarReconstruction(
            a_simplex, a_reconstruction, a_cutting_plane_index,
            a_moments_to_return);
  }
};

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CUT_SIMPLEX_DRIVERS_H_
