// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_GEOMETRY_TYPE_TRAITS_H_
#define SRC_GEOMETRY_GENERAL_GEOMETRY_TYPE_TRAITS_H_

#include "src/geometry/half_edge_structures/half_edge_polytope.h"
#include "src/geometry/polygons/divided_polygon.h"
#include "src/geometry/polygons/polygon.h"
#include "src/geometry/polygons/tri.h"
#include "src/geometry/polyhedrons/capped_dodecahedron.h"
#include "src/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_LLLL.h"
#include "src/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_LLLT.h"
#include "src/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_LTLT.h"
#include "src/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_LLTT.h"
#include "src/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_LTTT.h"
#include "src/geometry/polyhedrons/capped_dodecahedron_variations/capped_dodecahedron_TTTT.h"
#include "src/geometry/polyhedrons/capped_octahedron_variations/capped_octahedron_LLL.h"
#include "src/geometry/polyhedrons/capped_octahedron_variations/capped_octahedron_LLT.h"
#include "src/geometry/polyhedrons/capped_octahedron_variations/capped_octahedron_LTT.h"
#include "src/geometry/polyhedrons/capped_octahedron_variations/capped_octahedron_TTT.h"
#include "src/geometry/polyhedrons/symmetric_decompositions/symmetric_tet.h"
#include "src/geometry/polyhedrons/symmetric_decompositions/symmetric_pyramid.h"
#include "src/geometry/polyhedrons/symmetric_decompositions/symmetric_triangular_prism.h"
#include "src/geometry/polyhedrons/symmetric_decompositions/symmetric_hexahedron.h"
#include "src/geometry/polyhedrons/general_polyhedron.h"
#include "src/geometry/polyhedrons/concave_box.h"
#include "src/geometry/polyhedrons/dodecahedron.h"
#include "src/geometry/polyhedrons/hexahedron.h"
#include "src/geometry/polyhedrons/polyhedron_24.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/geometry/polyhedrons/octahedron.h"
#include "src/geometry/polyhedrons/triangular_prism.h"
#include "src/geometry/polyhedrons/pyramid.h"
#include "src/geometry/polyhedrons/tet.h"
#include "src/parameters/defined_types.h"

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

#endif  // SRC_GEOMETRY_GENERAL_GEOMETRY_TYPE_TRAITS_H_
