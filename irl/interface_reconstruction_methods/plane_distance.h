// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_PLANE_DISTANCE_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_PLANE_DISTANCE_H_

#include <algorithm>
#include <cmath>
#include <iostream>

#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/geometry/polyhedrons/tet.h"
#include "irl/helpers/helper.h"
#include "irl/interface_reconstruction_methods/iterative_distance_solver.h"
#include "irl/optimization/bisection.h"
#include "irl/optimization/secant.h"
#include "irl/parameters/constants.h"
#include "irl/planar_reconstruction/planar_separator.h"

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

#include "irl/interface_reconstruction_methods/plane_distance.tpp"

#endif // IRL_INTERFACE_RECONSTRUCTION_METHODS_PLANE_DISTANCE_H_
