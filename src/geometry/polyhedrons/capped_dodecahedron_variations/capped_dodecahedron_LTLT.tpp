// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_CAPPED_DODECAHEDRON_VARIATIONS_CAPPED_DODECAHEDRON_LTLT_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_CAPPED_DODECAHEDRON_VARIATIONS_CAPPED_DODECAHEDRON_LTLT_TPP_

#include <cassert>

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace capped_dodecahedron_LTLT_triangulation {
static constexpr UnsignedIndex_t datum_index = 5;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 8>
    face_triangle_decomposition{{{7, 8, 4},
                                 {7, 6, 8},
                                 {0, 1, 2},
                                 {0, 2, 3},
                                 {2, 6, 7},
                                 {2, 7, 3},
                                 {0, 3, 7},
                                 {0, 7, 4}}};
}  // namespace capped_dodecahedron_LTLT_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType> CappedDodecahedron_LTLTSpecialization<
    Derived, VertexType>::generateHalfEdgeVersion(void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void CappedDodecahedron_LTLTSpecialization<Derived, VertexType>::
    setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(40, 9, 13);

  for (UnsignedIndex_t v = 0; v < 9; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 40> ending_vertex_mapping{
      {8, 6, 5, 4, 8, 5, 8, 4, 7, 6, 8, 7, 1, 2, 3, 0, 5, 1, 0, 4,
       5, 0, 1, 5, 2, 5, 6, 2, 6, 7, 2, 7, 3, 2, 3, 7, 0, 7, 4, 0}};
  static constexpr std::array<UnsignedIndex_t, 40> previous_half_edge_mapping{
      {2,  0,  1,  5,  3,  4,  8,  6,  7,  11, 9,  10, 15, 12,
       13, 14, 18, 16, 17, 21, 19, 20, 24, 22, 23, 27, 25, 26,
       30, 28, 29, 33, 31, 32, 36, 34, 35, 39, 37, 38}};
  static constexpr std::array<UnsignedIndex_t, 40> next_half_edge_mapping{
      {1,  2,  0,  4,  5,  3,  7,  8,  6,  10, 11, 9,  13, 14,
       15, 12, 17, 18, 16, 20, 21, 19, 23, 24, 22, 26, 27, 25,
       29, 30, 28, 32, 33, 31, 35, 36, 34, 38, 39, 37}};
  static constexpr std::array<UnsignedIndex_t, 40> face_mapping{
      {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,  4,  4,  4,  4,  5,  5,  5,  6,
       6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12, 12}};
  static constexpr std::array<UnsignedIndex_t, 40> opposite_half_edge_mapping{
      {5,  10, 26, 20, 7,  0,  11, 4,  38, 29, 1,  6,  18, 22,
       33, 34, 21, 23, 12, 39, 3,  16, 13, 17, 25, 24, 2,  28,
       27, 9,  31, 30, 35, 14, 15, 32, 37, 36, 8,  19}};
  for (UnsignedIndex_t n = 0;
       n < static_cast<UnsignedIndex_t>(ending_vertex_mapping.size()); ++n) {
    HalfEdgeType& current_half_edge = a_half_edge_version->getHalfEdge(n);
    current_half_edge = HalfEdgeType(
        &a_half_edge_version->getVertex(ending_vertex_mapping[n]),
        &a_half_edge_version->getHalfEdge(previous_half_edge_mapping[n]),
        &a_half_edge_version->getHalfEdge(next_half_edge_mapping[n]),
        &a_half_edge_version->getFace(face_mapping[n]));
    current_half_edge.setOppositeHalfEdge(
        &a_half_edge_version->getHalfEdge(opposite_half_edge_mapping[n]));
    current_half_edge.getFace()->setStartingHalfEdge(&current_half_edge);
    current_half_edge.getVertex()->setHalfEdge(&current_half_edge);
  }
}

template <class Derived, class VertexType>
constexpr UnsignedIndex_t CappedDodecahedron_LTLTSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) {
  return static_cast<UnsignedIndex_t>(
      capped_dodecahedron_LTLT_triangulation::face_triangle_decomposition
          .size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
CappedDodecahedron_LTLTSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
  assert(a_tet < CappedDodecahedron_LTLTSpecialization::
                     getNumberOfSimplicesInDecomposition());
  return {
      capped_dodecahedron_LTLT_triangulation::face_triangle_decomposition[a_tet]
                                                                         [0],
      capped_dodecahedron_LTLT_triangulation::face_triangle_decomposition[a_tet]
                                                                         [1],
      capped_dodecahedron_LTLT_triangulation::face_triangle_decomposition[a_tet]
                                                                         [2],
      capped_dodecahedron_LTLT_triangulation::datum_index};
}

template <class Derived, class VertexType>
ProxyTet<Derived> CappedDodecahedron_LTLTSpecialization<Derived, VertexType>::
    getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

}  // namespace IRL
#endif  // SRC_GEOMETRY_POLYHEDRONS_CAPPED_DODECAHEDRON_VARIATIONS_CAPPED_DODECAHEDRON_LTLT_TPP_
