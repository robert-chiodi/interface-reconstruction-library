// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_INTERFACE_RECONSTRUCTION_METHODS_PLANE_DISTANCE_H_
#define SRC_INTERFACE_RECONSTRUCTION_METHODS_PLANE_DISTANCE_H_

#include <algorithm>
#include <cmath>
#include <iostream>

#include "src/data_structures/stack_vector.h"
#include "src/generic_cutting/generic_cutting.h"
#include "src/geometry/general/normal.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/geometry/polyhedrons/tet.h"
#include "src/helpers/helper.h"
#include "src/interface_reconstruction_methods/iterative_distance_solver.h"
#include "src/optimization/bisection.h"
#include "src/optimization/secant.h"
#include "src/parameters/constants.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {
/// \file plane_distance.h
///
/// This file contains functions to calculate and set
/// the distance to the plane(s) in provided
/// reconstructions so that a given liquid volume fraction
/// is recreated on the given cell.

/// \brief Volume conserving distance-finding routine
/// for single-plane reconstructions.
///
/// This function uses the analytical calculation from
/// Scardovelli & Zaleski, JCP 164,228-247 (2000) to calculate
/// the distance to a single plane to result in a
/// reconstruction with a given volume fractiom. In this function,
/// the main worker is `getAlpha(...)`.
///
/// \param[in] a_rectangular_cuboid Cell for which the reconstruction
/// with a plane of `a_normal` will be for.
/// \param[in] a_volume_fraction Liquid volume fraction the
/// reconstruction should result in.
/// \param[in] a_normal Normal for a plane that needs the
/// reconstruction.
double findDistanceOnePlane(const RectangularCuboid& a_rectangular_cuboid,
                            const double a_volume_fraction,
                            const Normal& a_normal);

double findDistanceOnePlane(const Tet& a_tet, const double a_volume_fraction,
                            const Normal& a_normal);

template <class CellType, class ReconstructionType>
inline void runIterativeSolverForDistance(
    const CellType& a_cell, const double a_volume_fraction,
    ReconstructionType* a_reconstruction,
    const double a_volume_fraction_tolerance =
        global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);

template <class CellType, class VolumeFractionArrayType>
inline void runIterativeSolverForDistance(
    const CellType& a_cell, const VolumeFractionArrayType& a_volume_fraction,
    PlanarSeparatorPathGroup* a_reconstruction,
    const double a_volume_fraction_tolerance =
        global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);

template <class CellType, class ReconstructionType>
inline void runProgressiveDistanceSolver(
    const CellType& a_cell, const double a_volume_fraction,
    ReconstructionType* a_reconstruction,
    const double a_volume_fraction_tolerance =
        global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);

}  // namespace IRL

#include "src/interface_reconstruction_methods/plane_distance.tpp"

#endif  // SRC_INTERFACE_RECONSTRUCTION_METHODS_PLANE_DISTANCE_H_
