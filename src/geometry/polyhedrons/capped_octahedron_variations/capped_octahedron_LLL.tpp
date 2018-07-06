// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_CAPPED_OCTAHEDRON_VARIATIONS_CAPPED_OCTAHEDRON_LLL_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_CAPPED_OCTAHEDRON_VARIATIONS_CAPPED_OCTAHEDRON_LLL_TPP_

#include <cassert>

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace capped_octahedron_LLL_triangulation {
static constexpr UnsignedIndex_t datum_index = 3;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 5>
    face_triangle_decomposition{
        {{0, 1, 2}, {6, 5, 4}, {0, 4, 1}, {1, 4, 5}, {1, 5, 2}}};
}  // namespace capped_octahedron_LLL_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType> CappedOctahedron_LLLSpecialization<
    Derived, VertexType>::generateHalfEdgeVersion(void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void CappedOctahedron_LLLSpecialization<Derived, VertexType>::
    setHalfEdgeVersion(HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(30, 7, 10);

  for (UnsignedIndex_t v = 0; v < 7; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 30> ending_vertex_mapping{
      {1, 2, 0, 4, 3, 6, 3, 5, 6, 5, 4, 6, 4, 1, 0,
       3, 4, 0, 4, 5, 1, 5, 2, 1, 5, 3, 2, 3, 0, 2}};
  static constexpr std::array<UnsignedIndex_t, 30> previous_half_edge_mapping{
      {2,  0,  1,  5,  3,  4,  8,  6,  7,  11, 9,  10, 14, 12, 13,
       17, 15, 16, 20, 18, 19, 23, 21, 22, 26, 24, 25, 29, 27, 28}};
  static constexpr std::array<UnsignedIndex_t, 30> next_half_edge_mapping{
      {1,  2,  0,  4,  5,  3,  7,  8,  6,  10, 11, 9,  13, 14, 12,
       16, 17, 15, 19, 20, 18, 22, 23, 21, 25, 26, 24, 28, 29, 27}};
  static constexpr std::array<UnsignedIndex_t, 30> face_mapping{
      {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
       5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9}};
  static constexpr std::array<UnsignedIndex_t, 30> opposite_half_edge_mapping{
      {14, 23, 29, 11, 16, 6,  5,  25, 9, 8,  19, 3,  17, 18, 0,
       28, 4,  12, 13, 10, 21, 20, 24, 1, 22, 7,  27, 26, 15, 2}};
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
constexpr UnsignedIndex_t CappedOctahedron_LLLSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) {
  return static_cast<UnsignedIndex_t>(
      capped_octahedron_LLL_triangulation::face_triangle_decomposition.size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
CappedOctahedron_LLLSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
  assert(a_tet < CappedOctahedron_LLLSpecialization::
                     getNumberOfSimplicesInDecomposition());
  return {
      capped_octahedron_LLL_triangulation::face_triangle_decomposition[a_tet]
                                                                      [0],
      capped_octahedron_LLL_triangulation::face_triangle_decomposition[a_tet]
                                                                      [1],
      capped_octahedron_LLL_triangulation::face_triangle_decomposition[a_tet]
                                                                      [2],
      capped_octahedron_LLL_triangulation::datum_index};
}

template <class Derived, class VertexType>
ProxyTet<Derived> CappedOctahedron_LLLSpecialization<Derived, VertexType>::
    getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

}  // namespace IRL
#endif  // SRC_GEOMETRY_POLYHEDRONS_CAPPED_OCTAHEDRON_VARIATIONS_CAPPED_OCTAHEDRON_LLL_TPP_
