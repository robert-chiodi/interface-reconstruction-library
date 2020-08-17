// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_RECURSIVE_SIMPLEX_CUTTING_INITIALIZER_H_
#define SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_RECURSIVE_SIMPLEX_CUTTING_INITIALIZER_H_

/// \brief Decomposition of geometric objects into simplices that
/// are then intersected with planes in reconstructions to generate
/// groups of simplices that lay solely on one side of the plane.
///
/// This form of cutting relies on the ability to separate an object into
/// a collection of signed simplices. In the context of IRL, each object
/// should have the methods getNumberOfSimplicesInDecomposition(), which is
/// the number of signed simplices that are needed to represent the original
/// object, and getSimplexFromDecomposition(UnsignedIndex_t), which returns
/// one of the simplices. This simplex is then recursively divided by the
/// reconstruction into more simplices that lay above and below the planes in
/// the reconstruction, ultimately being left with simplices that are either
/// below or above the reconstruction. Once that is the case, the moments for
/// the simplices are then calculate and accumulated in the volume_moments
/// variable. What happens to the newly generated simplicies is dictated by the
/// reconstruction type and type of moments that need to be returned. These
/// different behaviors are created through the use of common static interfaces
/// for treating the simplices, mainly defined in cut_simplex.h, which
/// specializes the class in cut_simplex_drivers.h. Ultimately, there are two
/// processes that happen. The simplex is separated into two by a plane and new
/// simplices above and below it can be generated using the lookup tables in
/// lookup.h. The simplices underneath the plane in the reconstruction continue
/// to be divided until no planes are left, at which time it is treated
/// according to a function in handle_enclosed_volume.h. Volumes above the
/// planes are treated through functins in continue_dividing_volume.h, which
/// handles things such as distributing these volumes over cutting networks or
/// simply ignoring the volumes above the planes. The special cases of 2D
/// simplices (triangles) or 3D simplices (tetrahedra) are handled through
/// the SimplexWrapper<SimplexType>, which is a specialized struct given
/// unique behaviors to certain functions depending on the dimension. These
/// simplex wrappers are conveniently described in simplex_wrappers.h

#include "src/generic_cutting/recursive_simplex_cutting/cut_simplex.h"
#include "src/generic_cutting/recursive_simplex_cutting/simplex_wrapper.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class ReturnType, class EncompassingType, class ReconstructionType>
inline ReturnType cutThroughRecursiveSimplex(
    const EncompassingType& a_encompassing_polyhedron,
    const ReconstructionType& a_separating_reconstruction);

}  // namespace IRL

#include "src/generic_cutting/recursive_simplex_cutting/recursive_simplex_cutting_initializer.tpp"

#endif  // SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_RECURSIVE_SIMPLEX_CUTTING_INITIALIZER_H_
