// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_CUT_POLYGON_H_
#define IRL_GENERIC_CUTTING_CUT_POLYGON_H_

#include <algorithm>

#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/recursive_simplex_cutting/lookup_tables.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/geometry/polygons/tri.h"
#include "irl/geometry/polyhedrons/hexahedron.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/helpers/geometric_cutting_helpers.h"
#include "irl/parameters/defined_types.h"
#include "irl/planar_reconstruction/planar_separator.h"

namespace IRL {

/// \brief Return a convex Polygon that represents the plane in
/// `a_reconstruction[a_plane_to_make_polygon]` on `a_polyhedron`.
///
/// NOTE: This function does not set the centroid stored in
/// a DividedPolygon when PolygonType is DividedPolygon.
///
/// Template requirements for `PolygonType`:
/// - Must be either Polygon or DividedPolygon.
///
/// Template requirements for `PolyhedronType`:
/// - Must be a Polyhedron that inherits from BasePolyhedron. Must also be
///   be capable of generating a PlanarLocalizer representing the convex
///   polyhedron.
///
///
template <class PolygonType, class PolyhedronType>
PolygonType getPlanePolygonFromReconstruction(
    const PolyhedronType& a_polyhedron, const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon);

/// \brief Cut a supplied convex polyhedron by a given reconstruction
///  and return the surface area.
///
/// This function cuts `a_polyhedron` by `a_reconstruction` to obtain
/// the interfacial surface area in that cell. This makes heavy use
/// of large lookup tables held in `lookup_tables.h`.
///
/// NOTE: If `a_reconstruction` contains multiple planes, it will
/// also cut the planes by one another before calculating surface area.
///
/// Template requirements for `PolyhedronType`:
/// - Must be a Polyhedron that inherits from BasePolyhedron.
///
template <class PolyhedronType>
double getReconstructionSurfaceArea(const PolyhedronType& a_polyhedron,
                                    const PlanarSeparator& a_reconstruction);

}  // namespace IRL

#include "irl/generic_cutting/cut_polygon.tpp"

#endif // IRL_GENERIC_CUTTING_CUT_POLYGON_H_
