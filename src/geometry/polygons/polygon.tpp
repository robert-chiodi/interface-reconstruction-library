// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_POLYGON_TPP_
#define SRC_GEOMETRY_POLYGONS_POLYGON_TPP_

#include <float.h>
#include <cmath>
#include <numeric>

#include "src/geometry/polygons/tri.h"
#include "src/helpers/geometric_cutting_helpers.h"
#include "src/helpers/helper.h"
#include "src/helpers/mymath.h"

namespace IRL {

template <class Derived, class VertexType>
HalfEdgePolygon<typename PolygonSpecialization<Derived, VertexType>::pt_type>
PolygonSpecialization<Derived, VertexType>::generateHalfEdgeVersion(
    void) const {
  HalfEdgePolygon<pt_type> half_edge_version;
  setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolygonType>
void PolygonSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolygonType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolygonType::half_edge_type;
  using FaceType = typename HalfEdgePolygonType::face_type;

  a_half_edge_version->resize(2 * this->getNumberOfVertices(),
                              this->getNumberOfVertices(), 1);

  a_half_edge_version->setPlaneOfExistence(this->getPlaneOfExistence());

  for (UnsignedIndex_t v = 0; v < this->getNumberOfVertices(); ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  HalfEdgeType& first_half_edge = a_half_edge_version->getHalfEdge(0);
  first_half_edge = HalfEdgeType(
      &a_half_edge_version->getVertex(1),
      &a_half_edge_version->getHalfEdge(this->getNumberOfVertices() - 1),
      &a_half_edge_version->getHalfEdge(1), &a_half_edge_version->getFace(0));
  first_half_edge.getVertex()->setHalfEdge(&first_half_edge);

  for (UnsignedIndex_t n = 1; n < this->getNumberOfVertices(); ++n) {
    HalfEdgeType& current_half_edge = a_half_edge_version->getHalfEdge(n);
    current_half_edge = HalfEdgeType(
        &a_half_edge_version->getVertex((n + 1) % this->getNumberOfVertices()),
        &a_half_edge_version->getHalfEdge(n - 1),
        &a_half_edge_version->getHalfEdge((n + 1) %
                                          this->getNumberOfVertices()),
        &a_half_edge_version->getFace(0));
    current_half_edge.getVertex()->setHalfEdge(&current_half_edge);
  }
  a_half_edge_version->getFace(0).setStartingHalfEdge(
      &a_half_edge_version->getHalfEdge(0));

  // Now make external loop..
  HalfEdgeType& first_outer_half_edge =
      a_half_edge_version->getHalfEdge(this->getNumberOfVertices());
  first_outer_half_edge = HalfEdgeType(
      &a_half_edge_version->getVertex(0),
      &a_half_edge_version->getHalfEdge(2 * this->getNumberOfVertices() - 1),
      &a_half_edge_version->getHalfEdge(this->getNumberOfVertices() + 1),
      &getOpenBoundaryFace<FaceType>());
  first_outer_half_edge.getVertex()->setHalfEdge(&first_outer_half_edge);

  for (UnsignedIndex_t n = this->getNumberOfVertices() + 1;
       n < 2 * (this->getNumberOfVertices()); ++n) {
    HalfEdgeType& current_half_edge = a_half_edge_version->getHalfEdge(n);
    current_half_edge = HalfEdgeType(
        &a_half_edge_version->getVertex(2 * (this->getNumberOfVertices()) - n),
        &a_half_edge_version->getHalfEdge(n - 1),
        &a_half_edge_version->getHalfEdge(
            (n - this->getNumberOfVertices() + 1) %
                (this->getNumberOfVertices()) +
            this->getNumberOfVertices()),
        &getOpenBoundaryFace<FaceType>());
  }

  // Set opposites
  a_half_edge_version->getHalfEdge(0).setOppositeHalfEdge(
      &a_half_edge_version->getHalfEdge(this->getNumberOfVertices()));
  a_half_edge_version->getHalfEdge(this->getNumberOfVertices())
      .setOppositeHalfEdge(&a_half_edge_version->getHalfEdge(0));

  UnsignedIndex_t opposite_index_difference = 0;
  for (UnsignedIndex_t n = this->getNumberOfVertices() + 1;
       n < 2 * this->getNumberOfVertices(); ++n) {
    opposite_index_difference += 2;
    HalfEdgeType& outer_half_edge = a_half_edge_version->getHalfEdge(n);
    HalfEdgeType& inner_half_edge =
        a_half_edge_version->getHalfEdge(n - opposite_index_difference);
    outer_half_edge.setOppositeHalfEdge(&inner_half_edge);
    inner_half_edge.setOppositeHalfEdge(&outer_half_edge);
  }
}

template <class Derived, class VertexType>
UnsignedIndex_t
PolygonSpecialization<Derived, VertexType>::getNumberOfSimplicesInDecomposition(
    void) const {
  assert(this->getNumberOfVertices() == 0 || this->getNumberOfVertices() > 2);
  return this->getNumberOfVertices() != 0 ? this->getNumberOfVertices() - 2 : 0;
}

template <class Derived, class VertexType>
std::array<UnsignedIndex_t, 3>
PolygonSpecialization<Derived, VertexType>::getSimplexIndicesFromDecomposition(
    UnsignedIndex_t a_tri_number_to_get) const {
  assert(a_tri_number_to_get < this->getNumberOfSimplicesInDecomposition());
  return {0, a_tri_number_to_get + 1, a_tri_number_to_get + 2};
}

template <class Derived, class VertexType>
ProxyTri<Derived>
PolygonSpecialization<Derived, VertexType>::getSimplexFromDecomposition(
    const UnsignedIndex_t a_tri_number_to_get) const {
  assert(a_tri_number_to_get < this->getNumberOfSimplicesInDecomposition());
  return ProxyTri<Derived>(
      static_cast<const Derived&>(*this),
      this->getSimplexIndicesFromDecomposition(a_tri_number_to_get));
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYGONS_POLYGON_TPP_
