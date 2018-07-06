// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_RECURSIVE_SIMPLEX_CUTTING_INITIALIZER_TPP_
#define SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_RECURSIVE_SIMPLEX_CUTTING_INITIALIZER_TPP_

namespace IRL {

template <class SimplexType>
inline bool isSimplexSufficientlyLarge(const SimplexType& a_simplex);

//******************************************************************* //
//     Template function definitions placed below this.
//******************************************************************* //
template <class ReturnType, class EncompassingType, class ReconstructionType>
inline ReturnType cutThroughRecursiveSimplex(
    const EncompassingType& a_encompassing_polyhedron,
    const ReconstructionType& a_separating_reconstruction) {
  ReturnType volume_moments;
  for (UnsignedIndex_t n = 0;
       n < a_encompassing_polyhedron.getNumberOfSimplicesInDecomposition();
       ++n) {
    const auto simplex_to_decompose =
        a_encompassing_polyhedron.getSimplexFromDecomposition(n);
    // Start cutting at first plane, plane 0.
    if (isSimplexSufficientlyLarge(simplex_to_decompose)) {
      getVolumeMomentsForSimplex(simplex_to_decompose,
                                 a_separating_reconstruction, 0,
                                 &volume_moments);
    }
  }
  return volume_moments;
}

template <class SimplexType>
bool isSimplexSufficientlyLarge(const SimplexType& a_simplex) {
  return a_simplex.calculateAbsoluteVolume() >
         SimplexWrapper<SimplexType>::minimumAmountToTrack();
}

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_RECURSIVE_SIMPLEX_CUTTING_INITIALIZER_TPP_
