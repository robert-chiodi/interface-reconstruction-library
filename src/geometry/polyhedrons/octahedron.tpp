// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_OCTAHEDRON_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_OCTAHEDRON_TPP_

#include <cassert>

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace octahedron_triangulation {
static constexpr UnsignedIndex_t datum_index = 0;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 3>
    face_triangle_decomposition{{{4, 3, 5}, {4, 5, 2}, {4, 2, 1}}};
}  // namespace octahedron_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType>
OctahedronSpecialization<Derived, VertexType>::generateHalfEdgeVersion(
    void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void OctahedronSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(24, 6, 8);

  for (UnsignedIndex_t v = 0; v < 6; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 24> ending_vertex_mapping{
      {1, 2, 0, 3, 5, 4, 5, 2, 4, 2, 1, 4, 1, 0, 4, 0, 3, 4, 2, 5, 0, 5, 3, 0}};
  static constexpr std::array<UnsignedIndex_t, 24> previous_half_edge_mapping{
      {2,  0,  1,  5,  3,  4,  8,  6,  7,  11, 9,  10,
       14, 12, 13, 17, 15, 16, 20, 18, 19, 23, 21, 22}};
  static constexpr std::array<UnsignedIndex_t, 24> next_half_edge_mapping{
      {1,  2,  0,  4,  5,  3,  7,  8,  6,  10, 11, 9,
       13, 14, 12, 16, 17, 15, 19, 20, 18, 22, 23, 21}};
  static constexpr std::array<UnsignedIndex_t, 24> face_mapping{
      {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7}};
  static constexpr std::array<UnsignedIndex_t, 24> opposite_half_edge_mapping{
      {13, 10, 18, 17, 22, 6, 5, 19, 9,  8,  1, 12,
       11, 0,  15, 14, 23, 3, 2, 7,  21, 20, 4, 16}};
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
constexpr UnsignedIndex_t OctahedronSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) {
  return static_cast<UnsignedIndex_t>(
      octahedron_triangulation::face_triangle_decomposition.size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
OctahedronSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
  assert(a_tet <
         OctahedronSpecialization::getNumberOfSimplicesInDecomposition());
  return {octahedron_triangulation::face_triangle_decomposition[a_tet][0],
          octahedron_triangulation::face_triangle_decomposition[a_tet][1],
          octahedron_triangulation::face_triangle_decomposition[a_tet][2],
          octahedron_triangulation::datum_index};
}

template <class Derived, class VertexType>
ProxyTet<Derived>
OctahedronSpecialization<Derived, VertexType>::getSimplexFromDecomposition(
    const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

}  // namespace IRL
#endif  // SRC_GEOMETRY_POLYHEDRONS_OCTAHEDRON_TPP_
