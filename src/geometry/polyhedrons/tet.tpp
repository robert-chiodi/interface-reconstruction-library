// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_TET_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_TET_TPP_

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/polygons/tri.h"
#include "src/helpers/geometric_cutting_helpers.h"

namespace IRL {
template <class Derived, class VertexType>
constexpr UnsignedIndex_t
TetSpecialization<Derived, VertexType>::getNumberOfSimplicesInDecomposition(
    void) {
  return 1;
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
TetSpecialization<Derived, VertexType>::getSimplexIndicesFromDecomposition(
    const UnsignedIndex_t a_tet) {
  assert(a_tet < TetSpecialization::getNumberOfSimplicesInDecomposition());
  return {0, 1, 2, 3};
}

template <class Derived, class VertexType>
ProxyTet<Derived>
TetSpecialization<Derived, VertexType>::getSimplexFromDecomposition(
    const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType>
TetSpecialization<Derived, VertexType>::generateHalfEdgeVersion(void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
inline Volume TetSpecialization<Derived, VertexType>::calculateAbsoluteVolume(
    void) const {
  return std::fabs(this->calculateVolume());
}

template <class Derived, class VertexType>
inline Volume TetSpecialization<Derived, VertexType>::calculateVolume(
    void) const {
  Volume3D_Functor volume;
  volume((*this));
  return volume.getMoments();
}

template <class Derived, class VertexType>
inline double TetSpecialization<Derived, VertexType>::calculateSign(
    void) const {
  return std::copysign(1.0, this->calculateVolume());
}

template <class Derived, class VertexType>
inline Pt TetSpecialization<Derived, VertexType>::calculateCentroid(
    void) const {
  return 0.25 * Pt((*this)[0] + (*this)[1] + (*this)[2] + (*this)[3]);
}

template <class Derived, class VertexType>
inline VolumeMoments TetSpecialization<Derived, VertexType>::calculateMoments()
    const {
  VolumeMoments3D_Functor volume_moments;
  volume_moments((*this));
  return volume_moments.getMoments();
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void TetSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;
  a_half_edge_version->resize(12, 4, 4);

  for (UnsignedIndex_t v = 0; v < 4; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 12> ending_vertex_mapping{
      {1, 2, 0, 3, 1, 0, 2, 3, 0, 3, 2, 1}};
  static constexpr std::array<UnsignedIndex_t, 12> previous_half_edge_mapping{
      {2, 0, 1, 5, 3, 4, 8, 6, 7, 11, 9, 10}};
  static constexpr std::array<UnsignedIndex_t, 12> next_half_edge_mapping{
      {1, 2, 0, 4, 5, 3, 7, 8, 6, 10, 11, 9}};
  static constexpr std::array<UnsignedIndex_t, 12> face_mapping{
      {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3}};
  static constexpr std::array<UnsignedIndex_t, 12> opposite_half_edge_mapping{
      {5, 11, 6, 8, 9, 0, 2, 10, 3, 4, 7, 1}};
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

namespace tet_detail {
template <class TriProxyType>
static Plane planeFromTri(const TriProxyType& a_Tri,
                          const Pt& a_volume_centroid);

template <class TriProxyType>
Plane planeFromTri(const TriProxyType& a_Tri, const Pt& a_volume_centroid) {
  Normal face_normal = a_Tri.calculateNormal();
  makeNormalFaceOutwardsFromVolume(a_volume_centroid, a_Tri.calculateCentroid(),
                                   &face_normal);
  return Plane(face_normal, face_normal * a_Tri[0]);
}

}  // namespace tet_detail

template <class Derived, class VertexType>
PlanarLocalizer TetSpecialization<Derived, VertexType>::getLocalizer(
    void) const {
  PlanarLocalizer reconstruction_to_return;
  reconstruction_to_return.setNumberOfPlanes(4);
  Pt volume_centroid = this->calculateCentroid();
  for (UnsignedIndex_t v = 0; v < 4; ++v) {
    const ProxyTri<Derived> tet_face = ProxyTri<Derived>::fromNoExistencePlane(
        static_cast<const Derived&>(*this), {v, (v + 1) % 4, (v + 2) % 4});
    reconstruction_to_return[v] =
        tet_detail::planeFromTri(tet_face, volume_centroid);
  }
  return reconstruction_to_return;
}

template <class VertexType>
template <class GeometryType>
StoredTet<VertexType>::StoredTet(const ProxyTet<GeometryType>& a_tet_proxy)
    : StoredVertexAccess<StoredTet<VertexType>, VertexType, 4>{
          a_tet_proxy[0], a_tet_proxy[1], a_tet_proxy[2], a_tet_proxy[3]} {}

template <class VertexType>
template <class GeometryType>
StoredTet<VertexType>& StoredTet<VertexType>::operator=(
    const ProxyTet<GeometryType>& a_tet_proxy) {
  (*this) = StoredTet(a_tet_proxy);
  return (*this);
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYHEDRONS_TET_TPP_
