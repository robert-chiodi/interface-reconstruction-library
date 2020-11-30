// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_INTERFACE_RECONSTRUCTION_METHODS_VOLUME_FRACTION_MATCHING_H_
#define IRL_INTERFACE_RECONSTRUCTION_METHODS_VOLUME_FRACTION_MATCHING_H_

#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/helpers/helper.h"
#include "irl/interface_reconstruction_methods/plane_distance.h"
#include "irl/interface_reconstruction_methods/reconstruction_cleaning.h"
#include "irl/planar_reconstruction/planar_separator_path_group.h"
#include "irl/parameters/constants.h"

namespace IRL {

template <class CellType, class PlanarType>
inline void setDistanceToMatchVolumeFraction(
    const CellType& a_cell, const double a_volume_fraction,
    PlanarType* a_reconstruction,
    const double a_volume_fraction_tolerance =
        global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);

template <class CellType, class VolumeFractionArrayType>
inline void setGroupDistanceToMatchVolumeFraction(
    const CellType& a_cell, const VolumeFractionArrayType& a_volume_fraction,
    PlanarSeparatorPathGroup* a_reconstruction,
    const double a_volume_fraction_tolerance =
        global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);

/// \brief Sets distance in `a_reconstruction` to
/// result in given liquid volume fraction.
///
/// This function sets the distance in `a_reconstruction`
/// to lead to the liquid volume fraction `a_volume_fraction`
/// in the unit computational cell [-0.5,0.5]^3. The
/// distance that will be set will be in reference to
/// the center of the unit cell as 0.0. If
/// different coordinate system is needed, make sure to
/// correct after use. If a one-plane reconstruction,
/// the volume fraction will be recreated exactly through
/// an analytical calculation. If a two-plane reconstruction
/// the volume fraction will be recreated to within
/// `a_volume_fraction_tolerance` using Newton-Raphson
/// optimization and potentially bisection.
///
/// \param[in] a_reconstruction Pointer to PlanarSeparator object that
/// will have its plane's distances set.
/// \param[in] a_rectangular_cuboid Cell for which the reconstruction
/// with a plane of `a_normal` will be for. Currently only used
/// for single-plane reconstruction. Two plane reconstruction
/// still assumed to be on unit cube.
/// \param[in] a_volume_fraction Liquid volume fraction
/// that the updated distance will recreate.
/// \param[in] a_volume_fraction_tolerance Tolerance allowed
/// in recreating `a_volume_fraction`.
template <class CellType, class PlanarType>
inline void setDistanceToMatchVolumeFractionPartialFill(
    const CellType& a_cell, const double a_volume_fraction,
    PlanarType* a_reconstruction,
    const double a_volume_fraction_tolerance =
        global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);

/// \brief Specialization for RectangularCuboids that calls Analytical distance
/// finding if a single plane.
template <class PlanarType>
inline void setDistanceToMatchVolumeFractionPartialFill(
    const RectangularCuboid& a_cell, const double a_volume_fraction,
    PlanarType* a_reconstruction,
    const double a_volume_fraction_tolerance =
        global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);

/// \brief Specialization for Tet  that calls Analytical distance
/// finding if a single plane.
template <class PlanarType>
inline void setDistanceToMatchVolumeFractionPartialFill(
    const Tet& a_cell, const double a_volume_fraction,
    PlanarType* a_reconstruction,
    const double a_volume_fraction_tolerance =
        global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);  

template <class CellType, class VolumeFractionArrayType>
inline void setGroupDistanceToMatchVolumeFractionPartialFill(
    const CellType& a_cell, const VolumeFractionArrayType& a_volume_fraction,
    PlanarSeparatorPathGroup* a_reconstruction,
    const double a_volume_fraction_tolerance =
        global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);
}  // namespace IRL

#include "irl/interface_reconstruction_methods/volume_fraction_matching.tpp"

#endif // IRL_INTERFACE_RECONSTRUCTION_METHODS_VOLUME_FRACTION_MATCHING_H_
