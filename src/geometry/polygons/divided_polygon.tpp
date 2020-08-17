// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_DIVIDED_POLYGON_TPP_
#define SRC_GEOMETRY_POLYGONS_DIVIDED_POLYGON_TPP_

#include "src/geometry/general/moment_calculation_through_simplices.h"

namespace IRL {

template <class Derived, class VertexType>
UnsignedIndex_t DividedPolygonSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) const {
  return this->getNumberOfVertices();
}

template <class Derived, class VertexType>
std::array<UnsignedIndex_t, 3>
DividedPolygonSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(
        UnsignedIndex_t a_tri_number_to_get) const {
  assert(a_tri_number_to_get < this->getNumberOfSimplicesInDecomposition());
  ++a_tri_number_to_get;
  return a_tri_number_to_get != this->getNumberOfSimplicesInDecomposition()
             ? std::array<UnsignedIndex_t, 3>(
                   {0, a_tri_number_to_get, a_tri_number_to_get + 1})
             : std::array<UnsignedIndex_t, 3>({0, a_tri_number_to_get, 1});
}

template <class Derived, class VertexType>
ProxyTri<MaskStripper<Derived>>
DividedPolygonSpecialization<Derived, VertexType>::getSimplexFromDecomposition(
    UnsignedIndex_t a_tri_number_to_get) const {
  assert(a_tri_number_to_get < this->getNumberOfSimplicesInDecomposition());
  return ProxyTri<MaskStripper<Derived>>(
      static_cast<const MaskStripper<Derived>&>(*this),
      this->getSimplexIndicesFromDecomposition(a_tri_number_to_get));
}

template <class Derived, class VertexType>
HalfEdgePolygon<VertexType>
DividedPolygonSpecialization<Derived, VertexType>::generateHalfEdgeVersion(
    void) const {
  HalfEdgePolygon<VertexType> half_edge_version;
  setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolygonType>
void DividedPolygonSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolygonType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolygonType::half_edge_type;
  using FaceType = typename HalfEdgePolygonType::face_type;

  a_half_edge_version->resize(4 * (this->getNumberOfVertices()),
                              this->getNumberOfVertices() + 1,
                              this->getNumberOfVertices());
  a_half_edge_version->setPlaneOfExistence(this->getPlaneOfExistence());

  for (UnsignedIndex_t v = 0; v < this->getNumberOfVertices(); ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }
  a_half_edge_version->getVertex(this->getNumberOfVertices())
      .setLocation(this->calculateCentroid());

  // Make polygon shell first
  HalfEdgeType& first_half_edge = a_half_edge_version->getHalfEdge(0);
  first_half_edge = HalfEdgeType(
      &a_half_edge_version->getVertex(1),
      &a_half_edge_version->getHalfEdge(this->getNumberOfVertices() - 1),
      &a_half_edge_version->getHalfEdge(1), &a_half_edge_version->getFace(0));
  first_half_edge.getVertex()->setHalfEdge(&first_half_edge);
  a_half_edge_version->getFace(0).setStartingHalfEdge(&first_half_edge);

  for (UnsignedIndex_t n = 1; n < this->getNumberOfVertices(); ++n) {
    HalfEdgeType& current_half_edge = a_half_edge_version->getHalfEdge(n);
    current_half_edge = HalfEdgeType(
        &a_half_edge_version->getVertex((n + 1) % this->getNumberOfVertices()),
        &a_half_edge_version->getHalfEdge(n - 1),
        &a_half_edge_version->getHalfEdge((n + 1) %
                                          this->getNumberOfVertices()),
        &a_half_edge_version->getFace(n));
    current_half_edge.getVertex()->setHalfEdge(&current_half_edge);
    current_half_edge.getFace()->setStartingHalfEdge(&current_half_edge);
  }

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

  // Now make internal triangularization
  auto& centroid_vertex =
      a_half_edge_version->getVertex(this->getNumberOfVertices());

  UnsignedIndex_t first_new_half_edge_index = 2 * this->getNumberOfVertices();
  for (UnsignedIndex_t n = 0; n < this->getNumberOfVertices(); ++n) {
    UnsignedIndex_t second_new_half_edge_index = first_new_half_edge_index + 1;
    HalfEdgeType& first_new_half_edge =
        a_half_edge_version->getHalfEdge(first_new_half_edge_index);
    HalfEdgeType& second_new_half_edge =
        a_half_edge_version->getHalfEdge(second_new_half_edge_index);
    HalfEdgeType& shell_half_edge = a_half_edge_version->getHalfEdge(n);

    first_new_half_edge = HalfEdgeType(&centroid_vertex, nullptr, nullptr,
                                       &a_half_edge_version->getFace(n));

    second_new_half_edge =
        HalfEdgeType(shell_half_edge.getPreviousVertex(), nullptr, nullptr,
                     &a_half_edge_version->getFace(n));

    doubleLinkHalfEdges(&shell_half_edge, &first_new_half_edge);
    doubleLinkHalfEdges(&first_new_half_edge, &second_new_half_edge);
    doubleLinkHalfEdges(&second_new_half_edge, &shell_half_edge);
    centroid_vertex.setHalfEdge(&first_new_half_edge);

    first_new_half_edge_index += 2;
  }

  // Now set opposite for internal parts
  first_new_half_edge_index = 2 * this->getNumberOfVertices();
  a_half_edge_version->getHalfEdge(first_new_half_edge_index + 1)
      .setOppositeHalfEdge(&a_half_edge_version->getHalfEdge(
          4 * (this->getNumberOfVertices()) - 2));
  a_half_edge_version->getHalfEdge(4 * (this->getNumberOfVertices()) - 2)
      .setOppositeHalfEdge(
          &a_half_edge_version->getHalfEdge(first_new_half_edge_index + 1));
  for (UnsignedIndex_t f = 0; f < this->getNumberOfVertices() - 1; ++f) {
    HalfEdgeType& first_new_half_edge =
        a_half_edge_version->getHalfEdge(first_new_half_edge_index);
    HalfEdgeType& opposite_half_edge =
        a_half_edge_version->getHalfEdge(first_new_half_edge_index + 3);
    first_new_half_edge.setOppositeHalfEdge(&opposite_half_edge);
    opposite_half_edge.setOppositeHalfEdge(&first_new_half_edge);

    first_new_half_edge_index += 2;
  }
}

template <class Derived, class VertexType>
void DividedPolygonSpecialization<Derived, VertexType>::resetCentroid(void) {
  this->getCentroid() = this->calculatePolygonCentroid();
}

template <class Derived, class VertexType>
Pt DividedPolygonSpecialization<Derived, VertexType>::calculateCentroid(
    void) const {
  return this->getCentroid();
}

template <class Derived, class VertexType>
VolumeMoments
DividedPolygonSpecialization<Derived, VertexType>::calculateMoments(
    void) const {
  Volume volume = this->calculateVolume();
  return {volume, volume * this->calculateCentroid()};
}

template <class Derived, class VertexType>
VertexType& DividedPolygonSpecialization<Derived, VertexType>::getCentroid(
    void) {
  return this->getDerived().unmaskedAccess(0);
}
template <class Derived, class VertexType>
const VertexType&
DividedPolygonSpecialization<Derived, VertexType>::getCentroid(void) const {
  return this->getDerived().unmaskedAccess(0);
}

template <class Derived, class VertexType>
Pt DividedPolygonSpecialization<Derived, VertexType>::calculatePolygonCentroid(
    void) const {
  return GeneralPolygon<Derived, VertexType,
                        ProxyTri<MaskStripper<Derived>>>::calculateCentroid();
}

template <class VertexType>
template <class AnyDerivedClass, class AnySimplexType>
ExpandableDividedPolygon<VertexType>
ExpandableDividedPolygon<VertexType>::fromPolygon(
    const GeneralPolygon<AnyDerivedClass, VertexType, AnySimplexType>&
        a_polygon) {
  return ExpandableDividedPolygon(a_polygon);
}
template <class VertexType>
template <class AnyDerivedClass, class AnySimplexType>
ExpandableDividedPolygon<VertexType>::ExpandableDividedPolygon(
    const GeneralPolygon<AnyDerivedClass, VertexType, AnySimplexType>&
        a_polygon) {
  this->setNumberOfVertices(a_polygon.getNumberOfVertices());
  for (UnsignedIndex_t n = 0; n < a_polygon.getNumberOfVertices(); ++n) {
    (*this)[n] = a_polygon[n];
  }
  this->getCentroid() = a_polygon.calculateCentroid();
  this->setPlaneOfExistence(a_polygon.getPlaneOfExistence());
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYGONS_DIVIDED_POLYGON_TPP_
