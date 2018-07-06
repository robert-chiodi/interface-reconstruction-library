// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CUT_SIMPLEX_H_
#define SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CUT_SIMPLEX_H_

#include <type_traits>
#include <utility>

#include "src/generic_cutting/general/class_classifications.h"
#include "src/generic_cutting/recursive_simplex_cutting/continue_dividing_volume.h"
#include "src/generic_cutting/recursive_simplex_cutting/cut_simplex_drivers.h"
#include "src/generic_cutting/recursive_simplex_cutting/handle_enclosed_volume.h"
#include "src/helpers/SFINAE_boiler_plate.h"

namespace IRL {

template <class SimplexType, class ReconstructionType, class ReturnType>
void getVolumeMomentsForSimplex(const SimplexType& a_simplex,
                                const ReconstructionType& a_reconstruction,
                                const UnsignedIndex_t a_cutting_plane_index,
                                ReturnType* a_moments_to_return);

}  // namespace IRL

#include "src/generic_cutting/recursive_simplex_cutting/cut_simplex.tpp"

#endif  // SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CUT_SIMPLEX_H_
