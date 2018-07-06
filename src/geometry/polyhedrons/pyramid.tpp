// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_PYRAMID_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_PYRAMID_TPP_

#include <cassert>

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace pyramid_triangulation {
static constexpr UnsignedIndex_t datum_index = 4;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 2>
    face_triangle_decomposition{{{0, 1, 2}, {0, 2, 3}}};
}  // namespace pyramid_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType>
PyramidSpecialization<Derived, VertexType>::generateHalfEdgeVersion(
    void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void PyramidSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(16, 5, 5);

  for (UnsignedIndex_t v = 0; v < 5; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 16> ending_vertex_mapping{
      {1, 2, 3, 0, 3, 2, 4, 2, 1, 4, 1, 0, 4, 0, 3, 4}};
  static constexpr std::array<UnsignedIndex_t, 16> previous_half_edge_mapping{
      {3, 0, 1, 2, 6, 4, 5, 9, 7, 8, 12, 10, 11, 15, 13, 14}};
  static constexpr std::array<UnsignedIndex_t, 16> next_half_edge_mapping{
      {1, 2, 3, 0, 5, 6, 4, 8, 9, 7, 11, 12, 10, 14, 15, 13}};
  static constexpr std::array<UnsignedIndex_t, 16> face_mapping{
      {0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4}};
  static constexpr std::array<UnsignedIndex_t, 16> opposite_half_edge_mapping{
      {11, 8, 5, 14, 15, 2, 7, 6, 1, 10, 9, 0, 13, 12, 3, 4}};
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
constexpr UnsignedIndex_t
PyramidSpecialization<Derived, VertexType>::getNumberOfSimplicesInDecomposition(
    void) {
  return static_cast<UnsignedIndex_t>(
      pyramid_triangulation::face_triangle_decomposition.size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
PyramidSpecialization<Derived, VertexType>::getSimplexIndicesFromDecomposition(
    const UnsignedIndex_t a_tet) {
  assert(a_tet < PyramidSpecialization::getNumberOfSimplicesInDecomposition());
  return {pyramid_triangulation::face_triangle_decomposition[a_tet][0],
          pyramid_triangulation::face_triangle_decomposition[a_tet][1],
          pyramid_triangulation::face_triangle_decomposition[a_tet][2],
          pyramid_triangulation::datum_index};
}

template <class Derived, class VertexType>
ProxyTet<Derived>
PyramidSpecialization<Derived, VertexType>::getSimplexFromDecomposition(
    const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

}  // namespace IRL
#endif  // SRC_GEOMETRY_POLYHEDRONS_PYRAMID_TPP_
