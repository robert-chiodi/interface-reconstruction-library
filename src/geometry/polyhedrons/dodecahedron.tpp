// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_DODECAHEDRON_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_DODECAHEDRON_TPP_

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace dodecahedron_triangulation {
static constexpr UnsignedIndex_t datum_index = 5;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 6>
    face_triangle_decomposition{
        {{3, 0, 1}, {3, 1, 2}, {4, 3, 7}, {4, 0, 3}, {3, 6, 7}, {3, 2, 6}}};
}  // namespace dodecahedron_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType>
DodecahedronSpecialization<Derived, VertexType>::generateHalfEdgeVersion(
    void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void DodecahedronSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(36, 8, 12);

  for (UnsignedIndex_t v = 0; v < 8; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 36> ending_vertex_mapping{
      {7, 6, 5, 4, 7, 5, 0, 1, 3, 1, 2, 3, 3, 7, 4, 0, 3, 4,
       5, 6, 2, 1, 5, 2, 5, 1, 0, 4, 5, 0, 6, 7, 3, 2, 6, 3}};
  static constexpr std::array<UnsignedIndex_t, 36> previous_half_edge_mapping{
      {2,  0,  1,  5,  3,  4,  8,  6,  7,  11, 9,  10, 14, 12, 13, 17, 15, 16,
       20, 18, 19, 23, 21, 22, 26, 24, 25, 29, 27, 28, 32, 30, 31, 35, 33, 34}};
  static constexpr std::array<UnsignedIndex_t, 36> next_half_edge_mapping{
      {1,  2,  0,  4,  5,  3,  7,  8,  6,  10, 11, 9,  13, 14, 12, 16, 17, 15,
       19, 20, 18, 22, 23, 21, 25, 26, 24, 28, 29, 27, 31, 32, 30, 34, 35, 33}};
  static constexpr std::array<UnsignedIndex_t, 36> face_mapping{
      {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4,  4,  4,  5,  5,  5,
       6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11}};
  static constexpr std::array<UnsignedIndex_t, 36> opposite_half_edge_mapping{
      {5,  31, 19, 28, 14, 0,  16, 26, 9, 8,  21, 33, 17, 32, 4,  27, 6,  12,
       23, 2,  34, 10, 25, 18, 29, 22, 7, 15, 3,  24, 35, 1,  13, 11, 20, 30}};
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
constexpr UnsignedIndex_t DodecahedronSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) {
  return static_cast<UnsignedIndex_t>(
      dodecahedron_triangulation::face_triangle_decomposition.size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
DodecahedronSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
  assert(a_tet <
         DodecahedronSpecialization::getNumberOfSimplicesInDecomposition());
  return {dodecahedron_triangulation::face_triangle_decomposition[a_tet][0],
          dodecahedron_triangulation::face_triangle_decomposition[a_tet][1],
          dodecahedron_triangulation::face_triangle_decomposition[a_tet][2],
          dodecahedron_triangulation::datum_index};
}

template <class Derived, class VertexType>
ProxyTet<Derived>
DodecahedronSpecialization<Derived, VertexType>::getSimplexFromDecomposition(
    const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYHEDRONS_DODECAHEDRON_TPP_
