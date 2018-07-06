// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_SYMMETRIC_HEXAHEDRON_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_SYMMETRIC_HEXAHEDRON_TPP_

#include <cassert>

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/half_edge_structures/half_edge.h"

namespace IRL {

namespace symmetric_hexahedron_triangulation {
static constexpr UnsignedIndex_t datum_index = 0;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 18>
    face_triangle_decomposition{{{1, 2, 8},
                                 {2, 3, 8},
                                 {5, 1, 9},
                                 {4, 5, 9},
                                 {1, 5, 10},
                                 {5, 6, 10},
                                 {6, 2, 10},
                                 {2, 1, 10},
                                 {2, 6, 11},
                                 {6, 7, 11},
                                 {7, 3, 11},
                                 {3, 2, 11},
                                 {3, 7, 12},
                                 {7, 4, 12},
                                 {5, 4, 13},
                                 {4, 7, 13},
                                 {7, 6, 13},
                                 {6, 5, 13}}};
}  // namespace symmetric_hexahedron_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType>
SymmetricHexahedronSpecialization<Derived, VertexType>::generateHalfEdgeVersion(
    void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void SymmetricHexahedronSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(72, 14, 24);

  for (UnsignedIndex_t v = 0; v < 14; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 72> ending_vertex_mapping{
      {1, 8,  0, 2, 8,  1, 3, 8,  2, 0, 8,  3, 1, 9,  5, 0, 9,  1,
       4, 9,  0, 5, 9,  4, 5, 10, 1, 6, 10, 5, 2, 10, 6, 1, 10, 2,
       6, 11, 2, 7, 11, 6, 3, 11, 7, 2, 11, 3, 3, 12, 0, 7, 12, 3,
       4, 12, 7, 0, 12, 4, 4, 13, 5, 7, 13, 4, 6, 13, 7, 5, 13, 6}};
  static constexpr std::array<UnsignedIndex_t, 72> previous_half_edge_mapping{
      {2,  0,  1,  5,  3,  4,  8,  6,  7,  11, 9,  10, 14, 12, 13, 17, 15, 16,
       20, 18, 19, 23, 21, 22, 26, 24, 25, 29, 27, 28, 32, 30, 31, 35, 33, 34,
       38, 36, 37, 41, 39, 40, 44, 42, 43, 47, 45, 46, 50, 48, 49, 53, 51, 52,
       56, 54, 55, 59, 57, 58, 62, 60, 61, 65, 63, 64, 68, 66, 67, 71, 69, 70}};
  static constexpr std::array<UnsignedIndex_t, 72> next_half_edge_mapping{
      {1,  2,  0,  4,  5,  3,  7,  8,  6,  10, 11, 9,  13, 14, 12, 16, 17, 15,
       19, 20, 18, 22, 23, 21, 25, 26, 24, 28, 29, 27, 31, 32, 30, 34, 35, 33,
       37, 38, 36, 40, 41, 39, 43, 44, 42, 46, 47, 45, 49, 50, 48, 52, 53, 51,
       55, 56, 54, 58, 59, 57, 61, 62, 60, 64, 65, 63, 67, 68, 66, 70, 71, 69}};
  static constexpr std::array<UnsignedIndex_t, 72> face_mapping{
      {0,  0,  0,  1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  5,  5,  5,
       6,  6,  6,  7,  7,  7,  8,  8,  8,  9,  9,  9,  10, 10, 10, 11, 11, 11,
       12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 17,
       18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 22, 23, 23, 23}};
  static constexpr std::array<UnsignedIndex_t, 72> opposite_half_edge_mapping{
      {15, 5,  10, 33, 8,  1,  45, 11, 4,  48, 2,  7,  24, 17, 22, 0,  20, 13,
       57, 23, 16, 60, 14, 19, 12, 29, 34, 69, 32, 25, 36, 35, 28, 3,  26, 31,
       30, 41, 46, 66, 44, 37, 51, 47, 40, 6,  38, 43, 9,  53, 58, 42, 56, 49,
       63, 59, 52, 18, 50, 55, 21, 65, 70, 54, 68, 61, 39, 71, 64, 27, 62, 67}};
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
constexpr UnsignedIndex_t SymmetricHexahedronSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) {
  return static_cast<UnsignedIndex_t>(
      symmetric_hexahedron_triangulation::face_triangle_decomposition.size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
SymmetricHexahedronSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
  assert(
      a_tet <
      SymmetricHexahedronSpecialization::getNumberOfSimplicesInDecomposition());
  return {
      symmetric_hexahedron_triangulation::face_triangle_decomposition[a_tet][0],
      symmetric_hexahedron_triangulation::face_triangle_decomposition[a_tet][1],
      symmetric_hexahedron_triangulation::face_triangle_decomposition[a_tet][2],
      symmetric_hexahedron_triangulation::datum_index};
}

template <class Derived, class VertexType>
ProxyTet<Derived> SymmetricHexahedronSpecialization<Derived, VertexType>::
    getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

template <class Derived, class VertexType>
void SymmetricHexahedronSpecialization<Derived, VertexType>::adjustCapToMatchVolume(
    const Volume a_correct_volume) {
  Volume starting_volume = this->calculateVolume();
  Volume needed_change = a_correct_volume - starting_volume;

  for (UnsignedIndex_t v = 4; v < 8; ++v) {
    (*this)[v].getPt() -= (*this)[13].getPt();
  }

  // Four tets that will create additional volume are the two
  // indices below, vertex 13, and vertex 14.
  // Vertex 14 is set to be vertex 13 + connecting_line*t.
  // (And therefore vertex 14 is not explicitly needed)  
  // This allows volume of each tet to be written as
  // vol = 1.0/6.0 * connecting_line*t \cdot
  // (crossProduct(tet_base_vertices[n][0], tet_base_vertices[n][1]))
  static constexpr std::array<std::array<UnsignedIndex_t, 2>, 4>
      tet_base_vertices{{{5, 4}, {4, 7}, {7, 6}, {6, 5}}};

  // Correction where V_8 is moved to recreate supplied volume.
  // It is moved in the direction of V_8 - opposite_face_centroid.
  auto summed_cross_products = Pt::fromScalarConstant(0.0);
  for (UnsignedIndex_t n = 0; n < tet_base_vertices.size(); ++n) {
    summed_cross_products +=
        crossProduct((*this)[tet_base_vertices[n][0]].getPt(),
                     (*this)[tet_base_vertices[n][1]].getPt());
  }

  // Move Vertex[8] in direction perpendicular to the face.
  Normal connecting_line = Normal::fromPtNormalized(summed_cross_products);

  double adjustment_along_line =
      6.0 * needed_change / safelyTiny(magnitude(summed_cross_products));

  for (UnsignedIndex_t v = 4; v < 8; ++v) {
    (*this)[v].getPt() += (*this)[13].getPt();
  }
  (*this)[13].getPt() += Normal::toPt(connecting_line * adjustment_along_line);
  assert(std::fabs(this->calculateVolume() - a_correct_volume) < 1.0e-14);
}  

}  // namespace IRL
#endif  // SRC_GEOMETRY_POLYHEDRONS_SYMMETRIC_DECOMPOSITIONS_SYMMETRIC_HEXAHEDRON_TPP_
