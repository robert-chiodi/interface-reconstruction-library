// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_VOLUME_FRACTION_MATCHING_H_
#define SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_VOLUME_FRACTION_MATCHING_H_

#include "src/c_interface/geometry/polyhedrons/c_hexahedron.h"
#include "src/c_interface/geometry/polyhedrons/c_pyramid.h"
#include "src/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "src/c_interface/geometry/polyhedrons/c_tet.h"
#include "src/c_interface/geometry/polyhedrons/c_triangular_prism.h"
#include "src/c_interface/planar_reconstruction/c_separators.h"
#include "src/c_interface/planar_reconstruction/c_planar_separator_path_group.h"
#include "src/geometry/polyhedrons/hexahedron.h"
#include "src/geometry/polyhedrons/pyramid.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/geometry/polyhedrons/tet.h"
#include "src/geometry/polyhedrons/triangular_prism.h"
#include "src/interface_reconstruction_methods/volume_fraction_matching.h"
#include "src/planar_reconstruction/planar_separator.h"

extern "C" {

/// \brief Reset the distances to each plane in order to
/// recreate the `a_volume_fraction` for a provided RectangularCuboid.
///
/// This function is used to set the distance to each plane in
/// the PlanarSeparator in order to recreate the volume fraction
/// `a_volume_fraction` in a RectangularCuboid provided.
/// If it is a single plane PlanarSeparator in a RectangularCuboid, then the
/// analytical relation of Scardovelli & Zaleski (JCP, 2000) is used.
/// In all other cases, a non-linear root search is performed using the Secant
/// method, which is then passed to a bisection routine if an answer is not
/// initially found.
///
/// @param[in] a_cell Pointer to a RectangularCuboid
/// @param[in] a_volume_fraction Volume fraction that should
/// be matched once distance is set.
/// @param[in] a_reconstruction PlanarSeparator to set distances for.
///
void c_matchVolumeFraction_RectCub_PlanarSep_Default(
    const c_RectCub* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction);

/// \brief Same as
/// c_setDistanceToMatchVolumeFraction_RectangularCuboid_PlanarSeparator_DefaultTolerance
/// however with a given tolerance to satisfy when finding distance.
void c_matchVolumeFraction_RectCub_PlanarSep(
    const c_RectCub* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction, const double* a_volume_fraction_tolerance);

void c_matchVolumeFraction_Hex_PlanarSep_Default(
    const c_Hex* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction);

void c_matchVolumeFraction_Hex_PlanarSep(
    const c_Hex* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction, const double* a_volume_fraction_tolerance);

void c_matchVolumeFraction_Tet_PlanarSep_Default(
    const c_Tet* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction);

void c_matchVolumeFraction_Tet_PlanarSep(
    const c_Tet* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction, const double* a_volume_fraction_tolerance);

void c_matchVolumeFraction_TriPrism_PlanarSep_Default(
    const c_TriPrism* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction);

void c_matchVolumeFraction_TriPrism_PlanarSep(
    const c_TriPrism* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction, const double* a_volume_fraction_tolerance);

void c_matchVolumeFraction_Pyrmd_PlanarSep_Default(
    const c_Pyrmd* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction);

void c_matchVolumeFraction_Pyrmd_PlanarSep(
    const c_Pyrmd* a_cell, const double* a_volume_fraction,
    c_PlanarSep* a_reconstruction, const double* a_volume_fraction_tolerance);

void c_matchGroupVolumeFraction_Tet_PlanarSepPathGroup_Default(
    const c_Tet* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction);

void c_matchGroupVolumeFraction_Tet_PlanarSepPathGroup(
    const c_Tet* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction, const double* a_volume_fraction_tolerance);

void c_matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup_Default(
    const c_Pyrmd* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction);

void c_matchGroupVolumeFraction_Pyrmd_PlanarSepPathGroup(
    const c_Pyrmd* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction, const double* a_volume_fraction_tolerance);

void c_matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup_Default(
    const c_TriPrism* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction);

void c_matchGroupVolumeFraction_TriPrism_PlanarSepPathGroup(
    const c_TriPrism* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction, const double* a_volume_fraction_tolerance);

void c_matchGroupVolumeFraction_Hex_PlanarSepPathGroup_Default(
    const c_Hex* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction);

void c_matchGroupVolumeFraction_Hex_PlanarSepPathGroup(
    const c_Hex* a_cell, const int* a_size, const double* a_volume_fraction,
    c_PlanarSepPathGroup* a_reconstruction, const double* a_volume_fraction_tolerance);

}  // end extern C

#endif  // SRC_C_INTERFACE_INTERFACE_RECONSTRUCTION_METHODS_C_VOLUME_FRACTION_MATCHING_H_
