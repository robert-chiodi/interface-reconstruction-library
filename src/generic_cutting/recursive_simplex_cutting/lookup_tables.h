// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_LOOKUP_TABLES_H_
#define SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_LOOKUP_TABLES_H_

#include <limits>

#include "src/parameters/defined_types.h"

namespace IRL {

// To indicate element of table should not be used.
namespace {
constexpr LookupIndex_t NA = std::numeric_limits<LookupIndex_t>::max();
}

/// \file lookup_tables.h
///
/// This file contains lookup tables that are needed
/// to perform the geoemtric cutting operations.
/// They are used to determine resulting tets from cutting operations
/// and dividing cells into 5 tets. In most cases, the ordering of the numbers
/// in the lookup table below is VERY important and must be kept
/// self-consistent. It is best not to change them unless you are very
/// familiar with the routines and what each lookup table is doing and
/// how it is used.

namespace cut_tet_by_plane {

///  Tet vertex ordering, showing positive volume tet
///
///                    1                  y
///                  .                  |
///                / | \                |
///              /   |   \              |_____ x
///         3  /     |     \. 0        /
///           \      |    /           /
///             \    |  /              z
///                \ |/
///                    2
///
///
///

/// \brief The total number of resulting single-side tets after cutting a tet
/// by a plane.
///
/// This is the number of single-side tets that result from cutting a tet by a
/// plane. By single-side tet, we mean one that lays entirely below or above
/// the cutting plane. By doing this, the original tet is represented by a
/// collection of tets that lay entirely on a side of the plane. When this
/// plane represents a liquid-gas interface for instance, the new tets are now
/// fully in a single phase, allowing separate liquid and gas quantities to be
/// found for the originally phase-mixed tet.
///
/// The index is [cutting_case], which indicates the sign of each tet vertex.
static constexpr LookupIndex_t number_of_tets_after_cut[16] = {
    1, 4, 4, 6, 4, 6, 6, 4, 4, 6, 6, 4, 6, 4, 4, 1};
/// \brief Number of negative tets (d<=0) formed after cutting a tet by a
/// plane.
///
/// Using this array along with the total number of tets
/// (`number_of_tets_after_cut[]`) the positive tets can also be cycled
/// through. The first index for positive tet is tet_first_positive =
/// number_of_negative_tets_after_cut[case], where the total number of
/// positive tets that exists is number_of_positive_tets_after_cut[case] =
/// number_of_tets_after_cut[case] - number_of_negative_tets_after_cut[case].
///
/// The index is [cutting_case], which indicates the sign of each tet vertex.
static constexpr LookupIndex_t number_of_negative_tets_after_cut[16] = {
    1, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 1, 3, 1, 1, 0};

/// \brief The number of intersections between the plane and tet edges.
///
/// This is the number of intersections between the cutting plane and tet
/// edges. When a plane intersects a tet edge, a new vertex will be placed at
/// this location. The new vertices are then used when constructing new tets
/// that lay completely on one side of the cutting plane.
///
/// The index is [cutting_case], which indicates the sign of each tet vertex.
static constexpr LookupIndex_t number_of_new_vertices_after_cut[16] = {
    0, 3, 3, 4, 3, 4, 4, 3, 3, 4, 4, 3, 4, 3, 3, 0};

/// \brief The vertices used in interpolation along the tet edge to find the
/// point of intersection with a plane.
///
/// This is a list of vertices that are used to interpolate to the
/// intersection point between tet edges and planes to create the new vertices
/// referred to in `number_of_new_vertices_after_cut[]`.
///
/// The three indices from left to right are:
/// - [cutting_case]
/// - [vertex_on_the_edge_to_get]
/// - [tet_edge_that_has_an_intersection]
static constexpr LookupIndex_t cut_vertices[16][2][4] = {
    {{NA, NA, NA, NA}, {NA, NA, NA, NA}}, {{0, 0, 0, NA}, {1, 2, 3, NA}},
    {{1, 1, 1, NA}, {0, 2, 3, NA}},       {{0, 0, 1, 1}, {2, 3, 2, 3}},
    {{2, 2, 2, NA}, {0, 1, 3, NA}},       {{0, 0, 2, 2}, {1, 3, 1, 3}},
    {{1, 1, 2, 2}, {0, 3, 0, 3}},         {{0, 1, 2, NA}, {3, 3, 3, NA}},

    {{0, 1, 2, NA}, {3, 3, 3, NA}},       {{1, 1, 2, 2}, {0, 3, 0, 3}},
    {{0, 0, 2, 2}, {1, 3, 1, 3}},         {{2, 2, 2, NA}, {0, 1, 3, NA}},
    {{0, 0, 1, 1}, {2, 3, 2, 3}},         {{1, 1, 1, NA}, {0, 2, 3, NA}},
    {{0, 0, 0, NA}, {1, 2, 3, NA}},       {{NA, NA, NA, NA}, {NA, NA, NA, NA}}};

/// \brief The vertices used to create the new tets in
/// `number_of_tets_after_cut[]`
///
/// This is the vertices of the tet used to create new, single-side tets after
/// cutting by a plane. Here, there are 8 possible vertex numbers (0-7), due
/// to the 4 original tet vertices and the possibility of 4 additional
/// vertices being created during intersection with the plane.
///
/// The three indices from left to right are:
/// - [cutting_case]
/// - [tet_number_to_get]
/// - [vertices_that_make_the_tet]
static constexpr LookupIndex_t verts_for_tets[16][6][4] = {{{0, 1, 2, 3},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA}},

                                                           {{4, 6, 3, 2},
                                                            {4, 3, 1, 2},
                                                            {4, 5, 6, 2},
                                                            {6, 5, 4, 0},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA}},

                                                           {{4, 0, 3, 2},
                                                            {4, 3, 6, 2},
                                                            {4, 6, 5, 2},
                                                            {5, 6, 4, 1},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA}},

                                                           {{6, 4, 5, 3},
                                                            {6, 5, 7, 3},
                                                            {6, 2, 4, 3},
                                                            {5, 4, 6, 1},
                                                            {5, 6, 7, 1},
                                                            {5, 0, 4, 1}},

                                                           {{5, 4, 0, 3},
                                                            {5, 0, 1, 3},
                                                            {5, 6, 4, 3},
                                                            {4, 6, 5, 2},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA}},

                                                           {{6, 7, 5, 3},
                                                            {6, 5, 4, 3},
                                                            {6, 4, 1, 3},
                                                            {4, 5, 7, 2},
                                                            {4, 7, 6, 2},
                                                            {4, 0, 5, 2}},

                                                           {{4, 5, 7, 3},
                                                            {4, 7, 6, 3},
                                                            {4, 6, 0, 3},
                                                            {7, 5, 4, 1},
                                                            {7, 4, 6, 1},
                                                            {7, 6, 2, 1}},

                                                           {{4, 5, 6, 3},
                                                            {6, 4, 0, 1},
                                                            {6, 0, 2, 1},
                                                            {6, 5, 4, 1},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA}},

                                                           {{6, 4, 0, 1},
                                                            {6, 0, 2, 1},
                                                            {6, 5, 4, 1},
                                                            {4, 5, 6, 3},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA}},

                                                           {{7, 5, 4, 1},
                                                            {7, 4, 6, 1},
                                                            {7, 6, 2, 1},
                                                            {4, 5, 7, 3},
                                                            {4, 7, 6, 3},
                                                            {4, 6, 0, 3}},

                                                           {{4, 5, 7, 2},
                                                            {4, 7, 6, 2},
                                                            {4, 0, 5, 2},
                                                            {6, 7, 5, 3},
                                                            {6, 5, 4, 3},
                                                            {6, 4, 1, 3}},

                                                           {{4, 6, 5, 2},
                                                            {5, 4, 0, 3},
                                                            {5, 0, 1, 3},
                                                            {5, 6, 4, 3},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA}},

                                                           {{5, 4, 6, 1},
                                                            {5, 6, 7, 1},
                                                            {5, 0, 4, 1},
                                                            {6, 4, 5, 3},
                                                            {6, 5, 7, 3},
                                                            {6, 2, 4, 3}},

                                                           {{5, 6, 4, 1},
                                                            {4, 0, 3, 2},
                                                            {4, 3, 6, 2},
                                                            {4, 6, 5, 2},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA}},

                                                           {{6, 5, 4, 0},
                                                            {4, 6, 3, 2},
                                                            {4, 3, 1, 2},
                                                            {4, 5, 6, 2},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA}},

                                                           {{0, 1, 2, 3},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA},
                                                            {NA, NA, NA, NA}}};
}  // namespace cut_tet_by_plane

namespace cut_tri_by_plane {
/// \brief The total number of resulting single-side tris after cutting a tri
/// by a plane.
///
/// This is the number of single-side tris that result from cutting a tri by a
/// plane. By single-side tri, we mean one that lays entirely below or above
/// the cutting plane. By doing this, the original tri is represented by a
/// collection of tets that lay entirely on a side of the plane. When this
/// plane represents a liquid-gas interface for instance, the new tri are now
/// fully in a single phase, allowing separate liquid and gas quantities to be
/// found for the originally phase-mixed tri.
///
/// The index is [cutting_case], which indicates the sign of each tri vertex.
static constexpr LookupIndex_t number_of_tris_after_cut[8] = {1, 3, 3, 3,
                                                              3, 3, 3, 1};
/// \brief Number of negative tris (d<=0) formed after cutting a tet by a
/// plane.
///
/// Using this array along with the total number of tris
/// (`number_of_tris_after_cut[]`) the positive tris can also be cycled
/// through. The first index for positive tri is tri_first_positive =
/// number_of_negative_tris_after_cut[case], where the total number of
/// positive tris that exists is number_of_positive_tris_after_cut[case] =
/// number_of_tris_after_cut[case] - number_of_negative_tris_after_cut[case].
///
/// The index is [cutting_case], which indicates the sign of each tet vertex.
static constexpr LookupIndex_t number_of_negative_tris_after_cut[8] = {
    1, 2, 2, 1, 2, 1, 1, 0};

/// \brief The number of intersections between the plane and tri edges.
///
/// This is the number of intersections between the cutting plane and tet
/// edges. When a plane intersects a tri edge, a new vertex will be placed at
/// this location. The new vertices are then used when constructing new tris
/// that lay completely on one side of the cutting plane.
///
/// The index is [cutting_case], which indicates the sign of each tri vertex.
static constexpr LookupIndex_t number_of_new_vertices_after_cut[8] = {
    0, 2, 2, 2, 2, 2, 2, 0};

/// \brief The vertices used in interpolation along the tri edge to find the
/// point of intersection with a plane.
///
/// This is a list of vertices that are used to interpolate to the
/// intersection point between tri edges and planes to create the new vertices
/// referred to in `number_of_new_vertices_after_cut[]`.
///
/// The three indices from left to right are:
/// - [cutting_case]
/// - [vertex_on_the_edge_to_get]
/// - [tri_edge_that_has_an_intersection]
static constexpr LookupIndex_t cut_vertices[8][2][2] = {
    {{NA, NA}, {NA, NA}}, {{0, 0}, {1, 2}},    {{1, 1}, {0, 2}},
    {{0, 1}, {2, 2}},     {{0, 1}, {2, 2}},    {{1, 1}, {0, 2}},
    {{0, 0}, {1, 2}},     {{NA, NA}, {NA, NA}}};

/// \brief The vertices used to create the new tris in
/// `number_of_tris_after_cut[]`
///
/// This is the vertices of the tri used to create new, single-side tris after
/// cutting by a plane. Here, there are 5 possible vertex numbers (0-4), due
/// to the 3 original tri vertices and the possibility of 2 additional
/// vertices being created during intersection with the plane.
///
/// The three indices from left to right are:
/// - [cutting_case]
/// - [tri_number_to_get]
/// - [vertices_that_make_the_tri]
static constexpr LookupIndex_t verts_for_tris[8][3][3] = {
    {{0, 1, 2}, {NA, NA, NA}, {NA, NA, NA}},
    {{3, 1, 2}, {3, 2, 4}, {3, 4, 0}},
    {{3, 4, 2}, {3, 2, 0}, {3, 1, 4}},
    {{3, 4, 2}, {4, 3, 0}, {4, 0, 1}},
    {{4, 3, 0}, {4, 0, 1}, {3, 4, 2}},
    {{3, 1, 4}, {3, 4, 2}, {3, 2, 0}},
    {{3, 4, 0}, {3, 1, 2}, {3, 2, 4}},
    {{0, 1, 2}, {NA, NA, NA}, {NA, NA, NA}}};
}  // namespace cut_tri_by_plane

namespace cut_plane_by_cuboid {
/// \brief Number of intersections that occur for a plane
/// intersecting a rectangular cuboid, according to a unique case ID.
///
/// Point ordering on the rectangular cuboid is important here.
/// Point ordering is as follows. (x,y,z) of vertex
///
/// Vertex 0 : (+, -, -)
/// Vertex 1 : (+, +, -)
/// Vertex 2 : (+, +, +)
/// Vertex 3 : (+, -, +)
/// Vertex 4 : (-, -, -)
/// Vertex 5 : (-, +, -)
/// Vertex 6 : (-, +, +)
/// Vertex 7 : (-, -, +)
///
/// The index is [cutting_case], which indicates the sign of each tet vertex.
static constexpr LookupIndex_t number_of_plane_intersections_with_cuboid[256] =
    {0,  3,  3,  4,  3,  NA, 4,  5,  3,  4,  NA, 5,  4,  5,  5,  4,  3,  4,  NA,
     5,  NA, NA, NA, NA, NA, 5,  NA, 6,  NA, NA, NA, 5,  3,  NA, 4,  5,  NA, NA,
     5,  6,  NA, NA, NA, NA, NA, NA, NA, 5,  4,  5,  5,  4,  NA, NA, NA, 5,  NA,
     NA, NA, 5,  NA, NA, NA, 4,  3,  NA, NA, NA, 4,  NA, 5,  NA, NA, NA, NA, NA,
     5,  NA, 6,  5,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
     NA, 4,  NA, 5,  NA, 5,  NA, 4,  5,  NA, NA, NA, NA, NA, NA, 5,  4,  5,  NA,
     6,  5,  NA, NA, 5,  4,  NA, NA, NA, NA, NA, NA, NA, 3,  3,  NA, NA, NA, NA,
     NA, NA, NA, 4,  5,  NA, NA, 5,  6,  NA, 5,  4,  5,  NA, NA, NA, NA, NA, NA,
     5,  4,  NA, 5,  NA, 5,  NA, 4,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
     NA, NA, NA, NA, NA, 5,  6,  NA, 5,  NA, NA, NA, NA, NA, 5,  NA, 4,  NA, NA,
     NA, 3,  4,  NA, NA, NA, 5,  NA, NA, NA, 5,  NA, NA, NA, 4,  5,  5,  4,  5,
     NA, NA, NA, NA, NA, NA, NA, 6,  5,  NA, NA, 5,  4,  NA, 3,  5,  NA, NA, NA,
     6,  NA, 5,  NA, NA, NA, NA, NA, 5,  NA, 4,  3,  4,  5,  5,  4,  5,  NA, 4,
     3,  5,  4,  NA, 3,  4,  3,  3,  0};

/// \brief The vertices used in interpolation along the cuboid edge to find
/// the point of intersection with a plane.
///
/// This is a list of vertices that are used to interpolate to the
/// intersection point between edges of a cuboid and planes to create the new
/// vertices referred to in `number_of_plane_intersections_with_cuboid[]`.
///
/// The three indices from left to right are:
/// - [cutting_case]
/// - [vertex_on_the_edge_to_get]
/// - [cuboid_edge_vertex_that_has_the_intersection]
static constexpr LookupIndex_t cuboid_cut_vertices[256][2][6] = {
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 3, 0, NA, NA, NA}, {1, 0, 4, NA, NA, NA}},
    {{0, 1, 1, NA, NA, NA}, {1, 5, 2, NA, NA, NA}},
    {{1, 3, 0, 1, NA, NA}, {2, 0, 4, 5, NA, NA}},
    {{1, 2, 2, NA, NA, NA}, {2, 6, 3, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 1, 2, 2, NA, NA}, {1, 5, 6, 3, NA, NA}},
    {{2, 3, 0, 1, 2, NA}, {3, 0, 4, 5, 6, NA}},
    {{2, 3, 3, NA, NA, NA}, {3, 7, 0, NA, NA, NA}},
    {{0, 2, 3, 0, NA, NA}, {1, 3, 7, 4, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 2, 3, 0, 1, NA}, {2, 3, 7, 4, 5, NA}},
    {{1, 2, 3, 3, NA, NA}, {2, 6, 7, 0, NA, NA}},
    {{0, 1, 2, 3, 0, NA}, {1, 2, 6, 7, 4, NA}},
    {{0, 1, 2, 3, 3, NA}, {1, 5, 6, 7, 0, NA}},
    {{0, 1, 2, 3, NA, NA}, {4, 5, 6, 7, NA, NA}},
    {{4, 0, 7, NA, NA, NA}, {5, 4, 4, NA, NA, NA}},
    {{0, 3, 7, 4, NA, NA}, {1, 0, 4, 5, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 3, 7, 4, 1, NA}, {2, 0, 4, 5, 5, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 2, 3, 7, 4, NA}, {1, 3, 7, 4, 5, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 2, 3, 7, 4, 1}, {2, 3, 7, 4, 5, 5}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{4, 1, 2, 3, 7, NA}, {5, 5, 6, 7, 4, NA}},
    {{4, 5, 1, NA, NA, NA}, {5, 6, 5, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 4, 5, 1, NA, NA}, {1, 5, 6, 2, NA, NA}},
    {{1, 3, 0, 4, 5, NA}, {2, 0, 4, 5, 6, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 4, 5, 2, 2, NA}, {1, 5, 6, 6, 3, NA}},
    {{2, 3, 0, 4, 5, 2}, {3, 0, 4, 5, 6, 6}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{4, 5, 2, 3, 0, NA}, {5, 6, 6, 7, 4, NA}},
    {{5, 1, 0, 7, NA, NA}, {6, 5, 4, 4, NA, NA}},
    {{0, 3, 7, 5, 1, NA}, {1, 0, 4, 6, 5, NA}},
    {{0, 0, 7, 5, 1, NA}, {1, 4, 4, 6, 2, NA}},
    {{1, 3, 7, 5, NA, NA}, {2, 0, 4, 6, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{2, 3, 7, 5, 2, NA}, {3, 0, 4, 6, 6, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 2, 3, 7, 5, NA}, {2, 3, 7, 4, 6, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{5, 2, 3, 7, NA, NA}, {6, 6, 7, 4, NA, NA}},
    {{5, 6, 2, NA, NA, NA}, {6, 7, 6, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 5, 6, 2, NA, NA}, {2, 6, 7, 3, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 1, 5, 6, 2, NA}, {1, 5, 6, 7, 3, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 5, 6, 3, 3, NA}, {2, 6, 7, 7, 0, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 1, 5, 6, 3, 3}, {1, 5, 6, 7, 7, 0}},
    {{5, 6, 3, 0, 1, NA}, {6, 7, 7, 4, 5, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{4, 6, 2, 1, NA, NA}, {5, 7, 6, 5, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 4, 6, 2, 1, NA}, {1, 5, 7, 6, 2, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 1, 4, 6, 2, NA}, {2, 5, 5, 7, 3, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 4, 6, 2, NA, NA}, {1, 5, 7, 3, NA, NA}},
    {{2, 3, 0, 4, 6, NA}, {3, 0, 4, 5, 7, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 4, 6, 3, 3, NA}, {1, 5, 7, 7, 0, NA}},
    {{4, 6, 3, 0, NA, NA}, {5, 7, 7, 4, NA, NA}},
    {{6, 2, 1, 0, 7, NA}, {7, 6, 5, 4, 4, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 0, 7, 6, 2, 1}, {1, 4, 4, 7, 6, 2}},
    {{1, 3, 7, 6, 2, NA}, {2, 0, 4, 7, 6, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 0, 7, 6, 2, NA}, {1, 4, 4, 7, 3, NA}},
    {{2, 3, 7, 6, NA, NA}, {3, 0, 4, 7, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{6, 3, 7, NA, NA, NA}, {7, 7, 4, NA, NA, NA}},
    {{6, 7, 3, NA, NA, NA}, {7, 4, 7, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{2, 6, 7, 3, NA, NA}, {3, 7, 4, 0, NA, NA}},
    {{0, 2, 6, 7, 0, NA}, {1, 3, 7, 4, 4, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 2, 6, 7, 3, NA}, {2, 6, 7, 4, 0, NA}},
    {{0, 1, 2, 6, 7, 0}, {1, 2, 6, 7, 4, 4}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{6, 7, 0, 1, 2, NA}, {7, 4, 4, 5, 6, NA}},
    {{4, 0, 3, 6, NA, NA}, {5, 4, 7, 7, NA, NA}},
    {{0, 3, 3, 6, 4, NA}, {1, 0, 7, 7, 5, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{2, 6, 4, 0, 3, NA}, {3, 7, 5, 4, 0, NA}},
    {{0, 2, 6, 4, NA, NA}, {1, 3, 7, 5, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 2, 6, 4, 1, NA}, {2, 3, 7, 5, 5, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 1, 2, 6, 4, NA}, {1, 2, 6, 7, 5, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{4, 1, 2, 6, NA, NA}, {5, 5, 6, 7, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{5, 1, 0, 3, 6, NA}, {6, 5, 4, 7, 7, NA}},
    {{0, 3, 3, 6, 5, 1}, {1, 0, 7, 7, 6, 5}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 3, 3, 6, 5, NA}, {2, 0, 7, 7, 6, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 2, 6, 5, 1, NA}, {1, 3, 7, 6, 5, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 2, 6, 5, NA, NA}, {2, 3, 7, 6, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{5, 2, 6, NA, NA, NA}, {6, 6, 7, NA, NA, NA}},
    {{5, 7, 3, 2, NA, NA}, {6, 4, 7, 6, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 5, 7, 3, 2, NA}, {2, 6, 4, 7, 3, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{2, 2, 5, 7, 3, NA}, {3, 6, 6, 4, 0, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 5, 7, 3, NA, NA}, {2, 6, 4, 0, NA, NA}},
    {{0, 1, 5, 7, 0, NA}, {1, 2, 6, 4, 4, NA}},
    {{0, 1, 5, 7, 3, NA}, {1, 5, 6, 4, 0, NA}},
    {{5, 7, 0, 1, NA, NA}, {6, 4, 4, 5, NA, NA}},
    {{4, 0, 3, 2, 5, NA}, {5, 4, 7, 6, 6, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{2, 2, 5, 4, 0, 3}, {3, 6, 6, 5, 4, 0}},
    {{0, 2, 2, 5, 4, NA}, {1, 3, 6, 6, 5, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 5, 4, 0, 3, NA}, {2, 6, 5, 4, 0, NA}},
    {{0, 1, 5, 4, NA, NA}, {1, 2, 6, 5, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{4, 1, 5, NA, NA, NA}, {5, 5, 6, NA, NA, NA}},
    {{4, 7, 3, 2, 1, NA}, {5, 4, 7, 6, 5, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 1, 4, 7, 3, 2}, {2, 5, 5, 4, 7, 3}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 4, 7, 3, 2, NA}, {1, 5, 4, 7, 3, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 1, 4, 7, 3, NA}, {2, 5, 5, 4, 0, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 4, 7, 3, NA, NA}, {1, 5, 4, 0, NA, NA}},
    {{4, 7, 0, NA, NA, NA}, {5, 4, 4, NA, NA, NA}},
    {{0, 3, 2, 1, NA, NA}, {4, 7, 6, 5, NA, NA}},
    {{0, 3, 3, 2, 1, NA}, {1, 0, 7, 6, 5, NA}},
    {{0, 0, 3, 2, 1, NA}, {1, 4, 7, 6, 2, NA}},
    {{1, 3, 3, 2, NA, NA}, {2, 0, 7, 6, NA, NA}},
    {{1, 1, 0, 3, 2, NA}, {2, 5, 4, 7, 3, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{0, 0, 3, 2, NA, NA}, {1, 4, 7, 3, NA, NA}},
    {{2, 3, 3, NA, NA, NA}, {3, 0, 7, NA, NA, NA}},
    {{2, 2, 1, 0, 3, NA}, {3, 6, 5, 4, 0, NA}},
    {{0, 2, 2, 1, NA, NA}, {1, 3, 6, 5, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}},
    {{1, 2, 2, NA, NA, NA}, {2, 3, 6, NA, NA, NA}},
    {{1, 1, 0, 3, NA, NA}, {2, 5, 4, 0, NA, NA}},
    {{0, 1, 1, NA, NA, NA}, {1, 2, 5, NA, NA, NA}},
    {{0, 0, 3, NA, NA, NA}, {1, 4, 0, NA, NA, NA}},
    {{NA, NA, NA, NA, NA, NA}, {NA, NA, NA, NA, NA, NA}}};
}  // namespace cut_plane_by_cuboid

namespace cut_plane_by_tet {

static constexpr LookupIndex_t number_of_plane_intersections_with_tet[16] = {
    0, 3, 3, 4, 3, 4, 4, 3, 3, 4, 4, 3, 4, 3, 3, 0};

static constexpr LookupIndex_t tet_cut_vertices[16][2][4] = {
    {{NA, NA, NA, NA}, {NA, NA, NA, NA}}, {{0, 0, 0, NA}, {1, 2, 3, NA}},
    {{1, 1, 1, NA}, {0, 3, 2, NA}},       {{0, 0, 1, 1}, {2, 3, 3, 2}},
    {{2, 2, 2, NA}, {0, 1, 3, NA}},       {{1, 1, 2, 3}, {0, 2, 3, 0}},
    {{1, 1, 2, 2}, {0, 3, 3, 0}},         {{0, 1, 2, NA}, {3, 3, 3, NA}},
    {{0, 2, 1, NA}, {3, 3, 3, NA}},       {{1, 2, 2, 1}, {0, 0, 3, 3}},
    {{0, 3, 2, 1}, {1, 0, 3, 2}},         {{2, 2, 2, NA}, {0, 3, 1, NA}},
    {{0, 1, 1, 0}, {2, 2, 3, 3}},         {{1, 1, 1, NA}, {0, 2, 3, NA}},
    {{0, 0, 0, NA}, {1, 3, 2, NA}},       {{NA, NA, NA, NA}, {NA, NA, NA, NA}}};

}  // namespace cut_plane_by_tet
}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_LOOKUP_TABLES_H_
