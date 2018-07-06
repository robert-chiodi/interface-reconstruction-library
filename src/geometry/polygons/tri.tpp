// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_TRI_TPP_
#define SRC_GEOMETRY_POLYGONS_TRI_TPP_

#include "src/helpers/geometric_cutting_helpers.h"

namespace IRL {

template <class Derived, class VertexType>
constexpr UnsignedIndex_t
TriSpecialization<Derived, VertexType>::getNumberOfSimplicesInDecomposition(
    void) {
  return 1;
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 3>
TriSpecialization<Derived, VertexType>::getSimplexIndicesFromDecomposition(
    UnsignedIndex_t a_tri_number_to_get) {
  assert(a_tri_number_to_get <
         TriSpecialization::getNumberOfSimplicesInDecomposition());
  return {0, 1, 2};
}

template <class Derived, class VertexType>
ProxyTri<Derived>
TriSpecialization<Derived, VertexType>::getSimplexFromDecomposition(
    const UnsignedIndex_t a_tri_number_to_get) const {
  assert(a_tri_number_to_get < this->getNumberOfSimplicesInDecomposition());
  return ProxyTri<Derived>(
      static_cast<const Derived&>(*this),
      this->getSimplexIndicesFromDecomposition(a_tri_number_to_get));
}

template <class Derived, class VertexType>
HalfEdgePolygon<typename TriSpecialization<Derived, VertexType>::pt_type>
TriSpecialization<Derived, VertexType>::generateHalfEdgeVersion(void) const {
  HalfEdgePolygon<pt_type> half_edge_version;
  setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolygonType>
void TriSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolygonType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolygonType::half_edge_type;
  using FaceType = typename HalfEdgePolygonType::face_type;

  a_half_edge_version->resize(2 * TriSpecialization::getNumberOfVertices(),
                              TriSpecialization::getNumberOfVertices(), 1);

  a_half_edge_version->setPlaneOfExistence(this->getPlaneOfExistence());

  for (UnsignedIndex_t v = 0; v < TriSpecialization::getNumberOfVertices();
       ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  a_half_edge_version->getHalfEdge(0) = HalfEdgeType(
      &a_half_edge_version->getVertex(1), &a_half_edge_version->getHalfEdge(2),
      &a_half_edge_version->getHalfEdge(1), &a_half_edge_version->getFace(0));
  a_half_edge_version->getHalfEdge(0).getVertex()->setHalfEdge(
      &a_half_edge_version->getHalfEdge(0));

  a_half_edge_version->getHalfEdge(1) = HalfEdgeType(
      &a_half_edge_version->getVertex(2), &a_half_edge_version->getHalfEdge(0),
      &a_half_edge_version->getHalfEdge(2), &a_half_edge_version->getFace(0));
  a_half_edge_version->getHalfEdge(1).getVertex()->setHalfEdge(
      &a_half_edge_version->getHalfEdge(1));

  a_half_edge_version->getHalfEdge(2) = HalfEdgeType(
      &a_half_edge_version->getVertex(0), &a_half_edge_version->getHalfEdge(1),
      &a_half_edge_version->getHalfEdge(0), &a_half_edge_version->getFace(0));
  a_half_edge_version->getHalfEdge(2).getVertex()->setHalfEdge(
      &a_half_edge_version->getHalfEdge(2));

  a_half_edge_version->getFace(0).setStartingHalfEdge(
      &a_half_edge_version->getHalfEdge(0));

  // External half-edges
  a_half_edge_version->getHalfEdge(3) = HalfEdgeType(
      &a_half_edge_version->getVertex(0), &a_half_edge_version->getHalfEdge(5),
      &a_half_edge_version->getHalfEdge(4), &getOpenBoundaryFace<FaceType>());

  a_half_edge_version->getHalfEdge(4) = HalfEdgeType(
      &a_half_edge_version->getVertex(2), &a_half_edge_version->getHalfEdge(3),
      &a_half_edge_version->getHalfEdge(5), &getOpenBoundaryFace<FaceType>());

  a_half_edge_version->getHalfEdge(5) = HalfEdgeType(
      &a_half_edge_version->getVertex(1), &a_half_edge_version->getHalfEdge(4),
      &a_half_edge_version->getHalfEdge(3), &getOpenBoundaryFace<FaceType>());

  // Set opposites
  a_half_edge_version->getHalfEdge(0).setOppositeHalfEdge(
      &a_half_edge_version->getHalfEdge(
          TriSpecialization::getNumberOfVertices()));
  a_half_edge_version->getHalfEdge(TriSpecialization::getNumberOfVertices())
      .setOppositeHalfEdge(&a_half_edge_version->getHalfEdge(0));

  UnsignedIndex_t opposite_index_difference = 0;
  for (UnsignedIndex_t n = TriSpecialization::getNumberOfVertices() + 1;
       n < 2 * TriSpecialization::getNumberOfVertices(); ++n) {
    opposite_index_difference += 2;
    HalfEdgeType& outer_half_edge = a_half_edge_version->getHalfEdge(n);
    HalfEdgeType& inner_half_edge =
        a_half_edge_version->getHalfEdge(n - opposite_index_difference);
    outer_half_edge.setOppositeHalfEdge(&inner_half_edge);
    inner_half_edge.setOppositeHalfEdge(&outer_half_edge);
  }
}

template <class Derived, class VertexType>
void TriSpecialization<Derived, VertexType>::reversePtOrdering(void) {
  std::swap((*this)[1], (*this)[2]);
}

template <class Derived, class VertexType>
Pt TriSpecialization<Derived, VertexType>::calculateCentroid(void) const {
  return Pt(((*this)[0] + (*this)[1] + (*this)[2])) / 3.0;
}

template <class Derived, class VertexType>
Volume TriSpecialization<Derived, VertexType>::calculateVolume(void) const {
  return this->calculateSign() * this->calculateAbsoluteVolume();
}

template <class Derived, class VertexType>
VolumeMomentsAndNormal
TriSpecialization<Derived, VertexType>::calculateVolumeMomentsAndNormal(
    void) const {
  VolumeMoments volume_moments = this->calculateMoments();
  return VolumeMomentsAndNormal(
      volume_moments, volume_moments.volume() * this->calculateNormal());
}

template <class Derived, class VertexType>
VolumeMoments TriSpecialization<Derived, VertexType>::calculateMoments(
    void) const {
  double volume = this->calculateVolume();
  return VolumeMoments(volume, volume * this->calculateCentroid());
}

template <class Derived, class VertexType>
double TriSpecialization<Derived, VertexType>::calculateSign(void) const {
  assert(magnitude(this->calculateNormal()) < 1.0e-10 ||
         magnitude(this->getPlaneOfExistence().normal()) > 0.9999);
  assert(std::fabs(this->getPlaneOfExistence().normal() *
                   this->calculateNormal()) > 0.9999);
  return this->getPlaneOfExistence().normal() * this->calculateNormal();
}

template <class Derived, class VertexType>
Volume TriSpecialization<Derived, VertexType>::calculateAbsoluteVolume(
    void) const {
  return 0.5 * magnitude(crossProduct((*this)[1] - (*this)[0],
		  	  	  	  	  	  	  	  (*this)[2] - (*this)[0]));
}

template <class VertexType>
template <class GeometryType>
StoredTri<VertexType>::StoredTri(const ProxyTri<GeometryType>& a_tri_proxy)
    : StaticStorageAndStoredPlane<StoredTri<VertexType>, VertexType, 3>(
          {a_tri_proxy[0], a_tri_proxy[1], a_tri_proxy[2]},
          a_tri_proxy.getPlaneOfExistence()) {}

template <class VertexType>
template <class GeometryType>
StoredTri<VertexType>& StoredTri<VertexType>::operator=(
    const ProxyTri<GeometryType>& a_tri_proxy) {
  (*this) = StoredTri(a_tri_proxy);
  return (*this);
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYGONS_TRI_TPP_
