// This file is part of the Interface Reconstruction Library (IRL)
// a library for interface reconstruction and computational geometry operations
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_CAPPED_OCTAHEDRON_VARIATIONS_CAPPED_OCTAHEDRON_CORRECTION_BASE_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_CAPPED_OCTAHEDRON_VARIATIONS_CAPPED_OCTAHEDRON_CORRECTION_BASE_TPP_

#include <cassert>

#include "src/geometry/general/moment_calculation_through_simplices.h"

 namespace IRL {

 template<class Derived, class VertexType>
 Derived& CappedOctahedronCorrectionBase<Derived, VertexType>::getDerived(void) {
   return static_cast<Derived&>(*this);
 }
 template<class Derived, class VertexType>
 const Derived& CappedOctahedronCorrectionBase<Derived, VertexType>::getDerived(
     void) const {
   return static_cast<const Derived&>(*this);
 }

template<class Derived, class VertexType>
VertexType& CappedOctahedronCorrectionBase<Derived, VertexType>::operator[](
    const UnsignedIndex_t a_index) {
  return this->getDerived().access(a_index);
}

template<class Derived, class VertexType>
const VertexType& CappedOctahedronCorrectionBase<Derived, VertexType>::
operator[](const UnsignedIndex_t a_index) const {
  return this->getDerived().access(a_index);
}

template <class Derived, class VertexType>
void CappedOctahedronCorrectionBase<Derived, VertexType>::
    adjustCapToMatchVolume(const Volume a_correct_volume) {
  Volume starting_volume = this->calculateVolume();
  Volume needed_change = a_correct_volume - starting_volume;

  for (UnsignedIndex_t v = 3; v < 6; ++v) {
    (*this)[v].getPt() -= (*this)[6].getPt();
  }

  // Three tets that will create additional volume are the two
  // indices below, vertex 6, and vertex 7.
  // Vertex 7 is set to be vertex 6 + connecting_line*t.
  // (and hence doesn't actually exist)
  // This allows volume of each tet to be written as
  // vol = 1.0/6.0 * connecting_line*t \cdot
  // (crossProduct(tet_base_vertices[n][0], tet_base_vertices[n][1]))
  static constexpr std::array<std::array<UnsignedIndex_t, 2>, 3>
      tet_base_vertices{{{4,3}, {3, 5}, {5, 4}}};

  // Correction where V_6 is moved to recreate supplied volume.
  // It is moved in the average orthogonal direction of the back face.
  auto summed_cross_products = Pt::fromScalarConstant(0.0);
  for (UnsignedIndex_t n = 0; n < 3; ++n) {
    summed_cross_products +=
        crossProduct((*this)[tet_base_vertices[n][0]].getPt(),
                     (*this)[tet_base_vertices[n][1]].getPt());
  }

  // Move Vertex[6] in direction perpendicular to the face.
  Normal connecting_line = Normal::fromPtNormalized(summed_cross_products);

  double adjustment_along_line =
      6.0 * needed_change / safelyTiny(magnitude(summed_cross_products));

  for (UnsignedIndex_t v = 3; v < 6; ++v) {
    (*this)[v].getPt() += (*this)[6].getPt();
  }
  (*this)[6].getPt() += Normal::toPt(connecting_line * adjustment_along_line);
  assert(std::fabs(this->calculateVolume() - a_correct_volume) < 1.0e-14);
}


 } // namespace IRL
#endif //SRC_GEOMETRY_POLYHEDRONS_CAPPED_OCTAHEDRON_VARIATIONS_CAPPED_OCTAHEDRON_CORRECTION_BASE_TPP_
