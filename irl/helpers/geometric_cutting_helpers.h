// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_HELPERS_GEOMETRIC_CUTTING_HELPERS_H_
#define IRL_HELPERS_GEOMETRIC_CUTTING_HELPERS_H_

#include <algorithm>
#include <array>
#include <cstring>
#include <numeric>

#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/pt.h"

namespace IRL {

/// \brief Calculates the distance between two points.
inline double distanceBetweenPts(const Pt& a_pt_1, const Pt& a_pt_2);

/// \brief Calculates the squared distance between two points.
inline double squaredDistanceBetweenPts(const Pt& a_pt_1, const Pt& a_pt_2);

/// \brief Finds the node in a ReconstructionLink network that the point is
/// internal to.
template <class ReconstructionType>
inline UnsignedIndex_t locatePt(
    const Pt& a_pt, const ReconstructionType& a_reconstruction_link);

/// \brief Correct normal to point outwards from the object, i.e.
/// in the same general direction as `a_surface_centroid - a_volume_centroid`.
inline void makeNormalFaceOutwardsFromVolume(const Pt& a_volume_centroid,
                                             const Pt& a_surface_centroid,
                                             Normal* a_normal);

/// \brief Given an array of distances to an intersecting plane
/// and the length of that array, a unique integer ID is calculated
/// based on the sign and index of each element in a_distances.
/// This is often used to then reference unique cases in a lookup table.
/// The ID is calculated as \f$
/// \sum_{i=0}^{i=\mathtt{a\_number\_of\_distances}}
/// 2^i*(0.5+\mathtt{sign}(0.5,\mathtt{a\_distances}[i])) \f$.

/// \param[in] a_distances Array of signed distances for the vertices.
/// \param[in] a_number_of_distances Length of the array a_distances.
template <class DistanceListType>
inline LookupIndex_t getGeometricCaseId(const DistanceListType& a_distances);

/// \brief Checks to see if a plane intersects a given cell.
///
/// Template requirements for `GeometryType`:
/// - Overload `operator[]` : Method that supplies
/// access to the underlying vertices as `Pt`s.
/// - `int getNumberOfVertices(void)`: Method that
/// returns the number of points in the geometry
/// `static const int max_number_of_vertices`:
/// Maximum number of pts that can be in the geometry,
/// used for allocation of a double array.
///
/// \param[in] a_plane Plane that may or may not intersect
/// `a_geometry`
/// \param[in] a_geometry Geometry on which to check if plane exists.
template <class GeometryType>
inline bool isPlaneIntersectingCell(const Plane& a_plane,
                                    const GeometryType& a_geometry);

/// \brief Checks to see if point is internal to PlanarReconstruction
///
/// This function determines if the point `a_pt` is internal to the
/// planar reconstruction or not. In the case of a single plane
/// in `a_reconstruction`, this would mean if the point is underneath
/// the plane.
///
/// Template requirements for `ReconstructionType`:
/// - Must be inherited from `PlanarReconstruction`,
/// such as a PlanarLocalizer or PlanarSeparator.
/// - double flip(void) : Returns -1.0 if cutting is flipped
/// (concave object) or 1.0 if cutting is not flipped (convex).
///
/// @param[in] a_pt Point for which it will be determined if it is internal or
/// not to `a_reconstruction`.
/// @param[in] a_reconstruction Reconstruction forming the volume
/// `a_pt` may or may not be internal to.
template <class ReconstructionType>
inline bool isPtInternal(const Pt& a_pt,
                         const ReconstructionType& a_reconstruction);

/// \brief Loops through the vertices of the supplied geometry and computes a
/// signed distance from `a_cutting_plane` for each.
///
/// This function computes the signed distance from `a_cutting_plane` to the
/// vertices in `a_geometry` using the `signedDistanceToPoint` method of
/// `a_cutting_plane`. `flip_cut` is needed in order to flip the sign of the
/// distance and cut for the correct phase when the reconstruction is comprised
/// of multiple planes. A supplied array of double[] is modified in the
/// function.
///
/// Template requirements for `GeometryType`:
/// - Overload `operator[]` : Method that supplies
/// access to the underlying vertices as `Pt`s.
/// - `int getNumberOfVertices(void)`: Method that
/// returns the number of points in the geometry
///
/// \param[in] a_tet The tet needing signed distances for its vertices.
/// \param[in] a_cutting_plane The plane that is being intersected with the tet.
/// \param[in] flip_cut Whether the sign should be flipped (indicating cutting
/// to obtain quantities for above the plane instead of below).
/// \param[out] a_distance_to_vertices Array of
/// double[GeometryType::max_number_of_vertices] to place computed distance for
/// vertices.
template <class GeometryType, class VertexDistanceList>
inline void signedDistanceToVertices(
    const GeometryType& a_geometry, const Plane& a_cutting_planem,
    VertexDistanceList* a_distance_to_vertices);

/// \brief Test if the object passed to it lays entirely
/// on one side of the reconstruction.
///
/// This function returns whether the polyhedron lays entirely
/// on one side of the reconstruction. Return values are as followed:
/// - 0 : Returned if polyhedron lays entirely internal to the reconstruction
/// (under the plane).
/// - 1 : Returned if polyhedron lays entirely external to the reconstruction
/// (over the plane).
/// - 2 : Returned if the reconstruction intersects the polyhedron.
template <class EncompassingType, class ReconstructionType>
inline int checkPureObject(const EncompassingType& a_encompassing_polyhedron,
                           const ReconstructionType& a_reconstruction);

}  // namespace IRL

#include "irl/helpers/geometric_cutting_helpers.tpp"

#endif // IRL_HELPERS_GEOMETRIC_CUTTING_HELPERS_H_
