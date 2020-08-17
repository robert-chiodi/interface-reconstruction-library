// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_CONCAVE_BOX_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_CONCAVE_BOX_TPP_

#include <cassert>

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace concave_box_triangulation {
static constexpr UnsignedIndex_t datum_index = 0;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 15>
    face_triangle_decomposition{{{2, 3, 8},
                                 {2, 8, 9},
                                 {3, 4, 7},
                                 {3, 7, 8},
                                 {4, 5, 6},
                                 {4, 6, 7},
                                 {5, 12, 13},
                                 {5, 13, 6},
                                 {12, 11, 10},
                                 {12, 10, 13},
                                 {1, 11, 12},
                                 {1, 12, 5},
                                 {1, 5, 4},
                                 {1, 4, 3},
                                 {1, 3, 2}}};
}  // namespace concave_box_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType>
ConcaveBoxSpecialization<Derived, VertexType>::generateHalfEdgeVersion(
    void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void ConcaveBoxSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(42, 14, 9);

  for (UnsignedIndex_t v = 0; v < 14; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 42> ending_vertex_mapping{
      {1,  2,  9,  0,  3, 8, 9, 2, 4,  7,  8,  3,  5,  6,
       7,  4,  12, 13, 6, 5, 1, 0, 10, 11, 11, 10, 13, 12,
       11, 12, 5,  4,  3, 2, 1, 9, 8,  7,  6,  13, 10, 0}};
  static constexpr std::array<UnsignedIndex_t, 42> previous_half_edge_mapping{
      {3,  0,  1,  2,  7,  4,  5,  6,  11, 8,  9,  10, 15, 12,
       13, 14, 19, 16, 17, 18, 23, 20, 21, 22, 27, 24, 25, 26,
       34, 28, 29, 30, 31, 32, 33, 41, 35, 36, 37, 38, 39, 40}};
  static constexpr std::array<UnsignedIndex_t, 42> next_half_edge_mapping{
      {1,  2,  3,  0,  5,  6,  7,  4,  9,  10, 11, 8,  13, 14,
       15, 12, 17, 18, 19, 16, 21, 22, 23, 20, 25, 26, 27, 24,
       29, 30, 31, 32, 33, 34, 28, 36, 37, 38, 39, 40, 41, 35}};
  static constexpr std::array<UnsignedIndex_t, 42> face_mapping{
      {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5,
       5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8}};
  static constexpr std::array<UnsignedIndex_t, 42> opposite_half_edge_mapping{
      {21, 34, 7,  35, 33, 11, 36, 2, 32, 15, 37, 5,  31, 19,
       38, 9,  30, 27, 39, 13, 28, 0, 41, 25, 29, 23, 40, 17,
       20, 24, 16, 12, 8,  4,  1,  3, 6,  10, 14, 18, 26, 22}};
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
constexpr UnsignedIndex_t ConcaveBoxSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) {
  return static_cast<UnsignedIndex_t>(
      concave_box_triangulation::face_triangle_decomposition.size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
ConcaveBoxSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
  assert(a_tet <
         ConcaveBoxSpecialization::getNumberOfSimplicesInDecomposition());
  return {concave_box_triangulation::face_triangle_decomposition[a_tet][0],
          concave_box_triangulation::face_triangle_decomposition[a_tet][1],
          concave_box_triangulation::face_triangle_decomposition[a_tet][2],
          concave_box_triangulation::datum_index};
}

template <class Derived, class VertexType>
ProxyTet<Derived>
ConcaveBoxSpecialization<Derived, VertexType>::getSimplexFromDecomposition(
    const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYHEDRONS_CONCAVE_BOX_TPP_
