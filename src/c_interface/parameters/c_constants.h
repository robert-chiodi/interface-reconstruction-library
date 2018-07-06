// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_PARAMETERS_C_CONSTANTS_H_
#define SRC_C_INTERFACE_PARAMETERS_C_CONSTANTS_H_

#include "src/parameters/constants.h"

extern "C" {

/// \file c_constants.h
///
/// These C-style funcions are
/// mapped to functions available in src/constants.h.
///
/// This file deals with functions that set global functions
/// involved in the IRL library.
///
/// Individual documentation for each
/// function is given alongside the function.

/// \brief Set VF_LOW and VF_HIGH defined in src/constants.cpp.
///
/// This function sets bounds on Volume Fraction (VF), setting VF_LOW and
/// VF_HIGH for all future computations. These are used as threshold values
/// to terminate some computations, especially during calculations of
/// interface reconstructions. VF_HIGH will automatically be set as 1.0 -
/// a_VF_low in order to preserve symmetry.
///
/// @param[in] a_VF_low Value to set VF_LOW.
void c_setVFBounds(const double* a_VF_low);

/// \brief Set the volume fraction tolerance for iterative distance finding.
///
/// This function sets the default volume fraction tolerance to be used
/// when an iterative distance finding routine is used. It will
/// always be the minimum of a_tolerance and `VF_LOW`.
///
/// @param[in] a_tolerance Default volume fraction tolerance to use
/// during iterative distance finding.
void c_setVFTolerance_IterativeDistanceFinding(const double* a_tolerance);

/// \brief Function to set `MINIMUM_VOLUME_TO_TRACK` defined in
/// src/constants.cpp.
///
/// This function sets `MINIMUM_VOLUME_TO_TRACK` to the value
/// `a_minimum_volume_to_track`. `MINIMUM_VOLUME_TO_TRACK`
/// is primarily used as the terminating condition for
/// the numerical integration and subdivision of polyhedra, where
/// sub-volumes less than `MINIMUM_VOLUME_TO_TRACK`
/// will be ignored.
///
/// @param[in] a_minimum_volume_to_track Value to set
/// `MINIMUM_VOLUME_TO_TRACK` to.
void c_setMinimumVolToTrack(const double* a_minimum_volume_to_track);

/// \brief Function to set `MINIMUM_SURFACE_AREA_TO_TRACK` defined in
/// src/constants.cpp.
///
/// This function sets `MINIMUM_SURFACE_AREA_TO_TRACK` to the value
/// `a_minimum_surface_area_to_track`. `MINIMUM_SURFACE_AREA_TO_TRACK`
/// is primarily used as the terminating condition for
/// the numerical integration and subdivision of polygons, where
/// sub-areas less than `MINIMUM_SURFACE_AREA_TO_TRACK`
/// will be ignored.
///
/// @param[in] a_minimum_surface_area_to_track Value to set
/// `MINIMUM_SURFACE_AREA_TO_TRACK` to.
void c_setMinimumSAToTrack(const double* a_minimum_surface_area_to_track);
}

#endif  // SRC_C_INTERFACE_PARAMETERS_C_CONSTANTS_H_
