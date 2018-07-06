// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_GENERIC_CUTTING_C_CUT_POLYGON_H_
#define SRC_C_INTERFACE_GENERIC_CUTTING_C_CUT_POLYGON_H_

#include "src/c_interface/geometry/polygons/c_divided_polygon.h"
#include "src/c_interface/geometry/polygons/c_polygon.h"
#include "src/c_interface/geometry/polyhedrons/c_rectangular_cuboid.h"
#include "src/c_interface/geometry/polyhedrons/c_tet.h"
#include "src/c_interface/geometry/polyhedrons/c_hexahedron.h"
#include "src/c_interface/planar_reconstruction/c_separators.h"
#include "src/generic_cutting/cut_polygon.h"
#include "src/geometry/polygons/tri.h"

extern "C" {

/// \file c_cut_polygon.h
///
/// These C-style funcions are
/// mapped to functions available in src/cut_polygon.h.
///
/// This file maps to functions that deal with intersecting polygons
/// with planes and calculating surface area from Polygons. This mostly means
/// the creation of Polygons from intersections of Planes and Polyhedra, or
/// intersection of Polygons with Planes to generate new (truncated) Polygons.
///
/// Individual documentation for each
/// function is given alongside the function.

/// \brief Create a Polygon by truncating a Plane from a PlanarSeparator
///  by a RectangularCuboid.
///
/// This function intersects the `a_plane_index` Plane of `a_separator` with a
/// `a_rectangular_cuboid` in order to generate a Polygon. If the
/// PlanarSeparator consists of multiple planes, the Polygon object will also be
/// intersected with the other planes in the PlanarSeparator.
///
/// @param[in] a_rectangular_cuboid Pointer to RectangularCuboid object
/// that will be used to truncate the plane.
/// @param[in] a_separator Pointer to PlanarSeparator object that
/// the plane from which a Polygon is being created is taken.
/// @param[in] a_plane_index Index of plane in `a_separator`
/// that the Polygon will be created from.
/// @param[out] a_polygon Pointer to Polygon object where the
/// created Polygon will be stored.
void c_getPoly_RectCub_Poly(const c_RectCub* a_rectangular_cuboid,
                            const c_PlanarSep* a_separator,
                            const int* a_plane_index, c_Poly* a_polygon);

void c_getPoly_Tet_Poly(const c_Tet* a_tet,
                               const c_PlanarSep* a_separator,
                               const int* a_plane_index,
                               c_Poly* a_divided_polygon);

void c_getPoly_Hex_Poly(const c_Hex* a_hexahedron,
                               const c_PlanarSep* a_separator,
                               const int* a_plane_index,
                               c_Poly* a_divided_polygon);

/// \brief Create a DividedPolygon by truncating a Plane from a PlanarSeparator
///  by a RectangularCuboid.
///
/// This function intersects the `a_plane_index` Plane of `a_separator` with a
/// `a_rectangular_cuboid` in order to generate a DividedPolygon. If the
/// PlanarSeparator consists of multiple planes, the Polygon object will also be
/// intersected with the other planes in the PlanarSeparator.
/// The centroid for the DividedPolygon is also updated before the function
/// returns.
///
/// @param[in] a_rectangular_cuboid Pointer to RectangularCuboid object
/// that will be used to truncate the plane.
/// @param[in] a_separator Pointer to PlanarSeparator object that
/// the plane from which a DividedPolygon is being created is taken.
/// @param[in] a_plane_index Index of plane in `a_separator`
/// that the DividedPolygon will be created from.
/// @param[out] a_polygon Pointer to DividedPolygon object where the
/// created Polygon will be stored.
void c_getPoly_RectCub_DivPoly(const c_RectCub* a_rectangular_cuboid,
                               const c_PlanarSep* a_separator,
                               const int* a_plane_index,
                               c_DivPoly* a_divided_polygon);

/// \brief Creates the interface polygons for the planes
/// in `a_separator` that exist in `a_rectangular_cuboid`
/// and returns the total area of these polygons.
///
/// This function creates Polygon objects for each plane
/// in `a_separator` that exists solely in `a_rectangular_cuboid`.
/// The area of these Polygon objects is summed and returned
/// from the function. If `a_separator` contains multiple planes
/// the polygons will also be intersected and truncated by them prior to the
/// calculation of the surface area.
///
/// @param[in] a_rectangular_cuboid Pointer to RectangularCuboid
/// object that will be used to truncate the planes in `a_separator`.
/// @param[in] a_separator Pointer to PlanarSeparator object
/// from which the planes will be taken.
double c_getSA_RectCub(const c_RectCub* a_rectangular_cuboid,
                       const c_PlanarSep* a_separator);

}  // end extern C

#endif  // SRC_C_INTERFACE_GENERIC_CUTTING_C_CUT_POLYGON_H_
