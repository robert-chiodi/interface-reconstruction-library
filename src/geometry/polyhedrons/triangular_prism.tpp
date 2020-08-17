// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_TRIANGULAR_PRISM_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_TRIANGULAR_PRISM_TPP_

#include <cassert>

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace triangular_prism_triangulation {
static constexpr UnsignedIndex_t datum_index = 0;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 3>
    face_triangle_decomposition{{{1, 4, 5}, {1, 5, 2}, {3, 5, 4}}};
}  // namespace triangular_prism_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType>
TriangularPrismSpecialization<Derived, VertexType>::generateHalfEdgeVersion(
    void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void TriangularPrismSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(18, 6, 5);

  for (UnsignedIndex_t v = 0; v < 6; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 18> ending_vertex_mapping{
      {1, 2, 0, 3, 4, 1, 0, 2, 5, 3, 0, 4, 5, 2, 1, 5, 4, 3}};
  static constexpr std::array<UnsignedIndex_t, 18> previous_half_edge_mapping{
      {2, 0, 1, 6, 3, 4, 5, 10, 7, 8, 9, 14, 11, 12, 13, 17, 15, 16}};
  static constexpr std::array<UnsignedIndex_t, 18> next_half_edge_mapping{
      {1, 2, 0, 4, 5, 6, 3, 8, 9, 10, 7, 12, 13, 14, 11, 16, 17, 15}};
  static constexpr std::array<UnsignedIndex_t, 18> face_mapping{
      {0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4}};
  static constexpr std::array<UnsignedIndex_t, 18> opposite_half_edge_mapping{
      {6, 14, 7, 10, 17, 11, 0, 2, 13, 15, 3, 5, 16, 8, 1, 9, 12, 4}};
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
constexpr UnsignedIndex_t TriangularPrismSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) {
  return static_cast<UnsignedIndex_t>(
      triangular_prism_triangulation::face_triangle_decomposition.size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
TriangularPrismSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
  assert(a_tet <
         TriangularPrismSpecialization::getNumberOfSimplicesInDecomposition());
  return {triangular_prism_triangulation::face_triangle_decomposition[a_tet][0],
          triangular_prism_triangulation::face_triangle_decomposition[a_tet][1],
          triangular_prism_triangulation::face_triangle_decomposition[a_tet][2],
          triangular_prism_triangulation::datum_index};
}

template <class Derived, class VertexType>
ProxyTet<Derived>
TriangularPrismSpecialization<Derived, VertexType>::getSimplexFromDecomposition(
    const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

}  // namespace IRL
#endif  // SRC_GEOMETRY_POLYHEDRONS_TRIANGULAR_PRISM_TPP_
