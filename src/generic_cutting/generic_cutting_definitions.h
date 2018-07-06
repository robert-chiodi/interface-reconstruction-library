// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_GENERIC_CUTTING_DEFINITIONS_H_
#define SRC_GENERIC_CUTTING_GENERIC_CUTTING_DEFINITIONS_H_

#include "src/generic_cutting/default_cutting_method.h"

namespace IRL {

template <class ReturnType, class EncompassingType>
__attribute__((pure)) __attribute__((hot)) inline ReturnType getVolumeMoments(
    const EncompassingType& a_encompassing_polyhedron);

template <class ReturnType, class CuttingMethod = DefaultCuttingMethod,
          class EncompassingType, class ReconstructionType>
__attribute__((pure)) __attribute__((hot)) inline ReturnType getVolumeMoments(
    const EncompassingType& a_encompassing_polyhedron,
    const ReconstructionType& a_reconstruction);

template <class ReturnType, class CuttingMethod = DefaultCuttingMethod,
          class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType>
__attribute__((hot)) inline ReturnType getVolumeMoments(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction);

template <class ReturnType, class CuttingMethod = DefaultCuttingMethod,
          class EncompassingType, class ReconstructionType>
__attribute__((pure)) __attribute__((hot)) inline ReturnType
getNormalizedVolumeMoments(const EncompassingType& a_encompassing_polyhedron,
                           const ReconstructionType& a_reconstruction);

template <class CuttingMethod = DefaultCuttingMethod, class EncompassingType,
          class ReconstructionType>
__attribute__((pure)) __attribute__((hot)) inline double getVolumeFraction(
    const EncompassingType& a_encompassing_polyhedron,
    const ReconstructionType& a_reconstruction);

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_GENERIC_CUTTING_DEFINITIONS_H_
