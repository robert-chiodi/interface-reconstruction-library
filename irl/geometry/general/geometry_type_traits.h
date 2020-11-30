// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_GEOMETRY_TYPE_TRAITS_H_
#define IRL_GEOMETRY_GENERAL_GEOMETRY_TYPE_TRAITS_H_

#include "irl/geometry/half_edge_structures/half_edge_polytope.h"
#include "irl/geometry/polygons/divided_polygon.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/geometry/polygons/tri.h"
#include "irl/geometry/polyhedrons/capped_dodecahedron.h"
#include "irl/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_LLLL.h"
#include "irl/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_LLLT.h"
#include "irl/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_LTLT.h"
#include "irl/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_LLTT.h"
#include "irl/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_LTTT.h"
#include "irl/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_TTTT.h"
#include "irl/geometry/polyhedrons/capped_octahedron_variations/capped_octahedron_LLL.h"
#include "irl/geometry/polyhedrons/capped_octahedron_variations/capped_octahedron_LLT.h"
#include "irl/geometry/polyhedrons/capped_octahedron_variations/capped_octahedron_LTT.h"
#include "irl/geometry/polyhedrons/capped_octahedron_variations/capped_octahedron_TTT.h"
#include "irl/geometry/polyhedrons/symmetric_decompositions/symmetric_tet.h"
#include "irl/geometry/polyhedrons/symmetric_decompositions/symmetric_pyramid.h"
#include "irl/geometry/polyhedrons/symmetric_decompositions/symmetric_triangular_prism.h"
#include "irl/geometry/polyhedrons/symmetric_decompositions/symmetric_hexahedron.h"
#include "irl/geometry/polyhedrons/general_polyhedron.h"
#include "irl/geometry/polyhedrons/concave_box.h"
#include "irl/geometry/polyhedrons/dodecahedron.h"
#include "irl/geometry/polyhedrons/hexahedron.h"
#include "irl/geometry/polyhedrons/polyhedron_24.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/geometry/polyhedrons/octahedron.h"
#include "irl/geometry/polyhedrons/triangular_prism.h"
#include "irl/geometry/polyhedrons/pyramid.h"
#include "irl/geometry/polyhedrons/tet.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class C>
struct is_polyhedron : std::false_type {};

template <class C>
struct is_polyhedron<const C> : is_polyhedron<C> {};

template <class VertexType>
struct is_polyhedron<StoredCappedDodecahedron<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredCappedDodecahedron_LLLL<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredCappedDodecahedron_LLLT<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredCappedDodecahedron_LTLT<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredCappedDodecahedron_LLTT<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredCappedDodecahedron_LTTT<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredCappedDodecahedron_TTTT<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredCappedOctahedron_LLL<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredCappedOctahedron_LLT<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredCappedOctahedron_LTT<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredCappedOctahedron_TTT<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredSymmetricTet<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredSymmetricPyramid<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredSymmetricTriangularPrism<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredSymmetricHexahedron<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredDodecahedron<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredHexahedron<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredRectangularCuboid<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredTet<VertexType>> : std::true_type {};

template <class GeometryType>
struct is_polyhedron<ProxyTet<GeometryType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredConcaveBox<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredPolyhedron24<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredTriangularPrism<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredOctahedron<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polyhedron<StoredPyramid<VertexType>> : std::true_type {};

template <class VertexType>
class SegmentedDecomposedPolyhedron;

template <class VertexType>
struct is_polyhedron<SegmentedDecomposedPolyhedron<VertexType>>
    : std::true_type {};

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
struct is_polyhedron<
    SegmentedHalfEdgePolyhedron<FaceType, VertexType, kMaxFaces, kMaxVertices>>
    : std::true_type {};

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
struct is_polyhedron<
  StoredGeneralPolyhedron<VertexType, kStaticAllocSize>>
    : std::true_type {};


// General Polyhedron
template <class C>
struct is_general_polyhedron : std::false_type {};

template <class C>
struct is_general_polyhedron<const C> : is_polyhedron<C> {};

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
struct is_general_polyhedron<
  StoredGeneralPolyhedron<VertexType, kStaticAllocSize>>
    : std::true_type {};

// Tetrahedron
template <class C>
struct is_tet : std::false_type {};

template <class C>
struct is_tet<const C> : is_tet<C> {};

template <class VertexType>
struct is_tet<StoredTet<VertexType>> : std::true_type {};

template <class GeometryType>
struct is_tet<ProxyTet<GeometryType>> : std::true_type {};

// Polygons
template <class C>
struct is_polygon : std::false_type {};

template <class C>
struct is_polygon<const C> : is_polygon<C> {};

template <class VertexType>
struct is_polygon<ExpandablePolygon<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polygon<ExpandableDividedPolygon<VertexType>> : std::true_type {};

template <class VertexType>
struct is_polygon<StoredTri<VertexType>> : std::true_type {};

template <class GeometryType>
struct is_polygon<ProxyTri<GeometryType>> : std::true_type {};

template <class VertexType>
class SegmentedDecomposedPolygon;

template <class VertexType>
struct is_polygon<SegmentedDecomposedPolygon<VertexType>> : std::true_type {};

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
struct is_polygon<
    SegmentedHalfEdgePolygon<FaceType, VertexType, kMaxFaces, kMaxVertices>>
    : std::true_type {};

// Triangle
template <class C>
struct is_tri : std::false_type {};

template <class C>
struct is_tri<const C> : is_tri<C> {};

template <class VertexType>
struct is_tri<StoredTri<VertexType>> : std::true_type {};

template <class GeometryType>
struct is_tri<ProxyTri<GeometryType>> : std::true_type {};

}  // namespace IRL

#endif // IRL_GEOMETRY_GENERAL_GEOMETRY_TYPE_TRAITS_H_
