// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PARAMETERS_CONSTANTS_H_
#define SRC_PARAMETERS_CONSTANTS_H_

#include <float.h>
#include <cmath>

#include "src/parameters/defined_types.h"

namespace IRL {

/// \brief Set VF_LOW and VF_HIGH defined in src/constants.cpp.
///
/// This function sets VF_LOW for all future computations. VF_LOW
/// is used as a lower threshold to terminate some computations,
/// especially during calculations of interface reconstructions.
/// VF_HIGH will automatically be set as 1.0 - a_VF_low in
/// order to preserve symmetry.
///
/// @param[in] a_VF_low Value to set VF_LOW.
void setVolumeFractionBounds(const double a_VF_low);

/// \brief Set the volume fraction tolerance for iterative distance finding.
///
/// This function sets the volume fraction tolerance to be used
/// when an iterative distance finding routine is used. It will
/// always be the minimum of a_tolerance and `VF_LOW`.
///
/// @param[in] a_tolerance Volume fraction tolerance to use
/// during iterative distance finding.
void setVolumeFractionTolerance(const double a_tolerance);

/// \brief Function to set `MINIMUM_VOLUME_TO_TRACK`
void setMinimumVolumeToTrack(const double a_minimum_volume_to_track);

/// \brief Function to set `MINIMUM_SURFACE_AREA_TO_TRACK`
void setMinimumSurfaceAreaToTrack(const double a_minimum_surface_area_to_track);

namespace global_constants {

/// \file constants.h
///
/// This file contains constants related to
/// reconstructions and the computational geometry
/// routines.

/// \brief Volume fraction below which VOF will be set to 0.0.
extern double VF_LOW;

/// \brief Volume fraction above which VOF will be set to 1.0.
extern double VF_HIGH;

/// \brief An arbitrarily large distance for planes of single-phase cells.
///
/// An arbitrarily large distance for a plane of a single cell. Used so
/// that no matter the normal, the cell will be seen as single phase
/// once cut computationally. The only requirement on for
/// `ARBITRARILY_LARGE_DISTANCE` is that it is larger than
/// the size of the domain (i.e. > Lx + Ly + Lz) if global datum
/// reconstructions are used. If only local reconstruction are used,
/// it only needs to be large enough to place the place outside of
/// the unit cell.
static constexpr double ARBITRARILY_LARGE_DISTANCE = 0.5 * DBL_MAX;

/// \brief Default tolerance to achieve in matching volume fraction during two
/// plane reconstructions.
extern double TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE;

/// \brief Maximum number of planes that will below in a PlanarSeparator
constexpr UnsignedIndex_t MAX_PLANAR_SEPARATOR_PLANES = 3;

/// \brief Maximum number of planes that will below in a PlanarLoalizer
constexpr UnsignedIndex_t MAX_PLANAR_LOCALIZER_PLANES = 8;

/// \brief Max number of planes in a PlanarReconstruction based object.
constexpr UnsignedIndex_t MAX_PLANES_IN_OBJECT =
    MAX_PLANAR_SEPARATOR_PLANES > MAX_PLANAR_LOCALIZER_PLANES
        ? MAX_PLANAR_SEPARATOR_PLANES
        : MAX_PLANAR_LOCALIZER_PLANES;

/// \brief Value that determines if two vectors are the same through \f$n_1
/// \cdot n_2 > SAMEVEC \f$
static constexpr double SAME_VEC = 1.0 - 1.0e-4;

/// \brief Minimum volume we care about keeping track of
extern double MINIMUM_VOLUME_TO_TRACK;

/// \brief Minimum surface area we care about keeping track of
extern double MINIMUM_SURFACE_AREA_TO_TRACK;

}  // namespace global_constants

}  // namespace IRL

#endif  // SRC_PARAMETERS_CONSTANTS_H_
