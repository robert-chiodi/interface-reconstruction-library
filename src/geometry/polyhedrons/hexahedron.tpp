// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_HEXAHEDRON_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_HEXAHEDRON_TPP_

namespace IRL {

namespace hexahedron_triangulation {
static constexpr UnsignedIndex_t datum_index = 4;
static constexpr std::array<std::array<UnsignedIndex_t, 3>, 6>
    face_triangle_decomposition{
        {{0, 1, 2}, {0, 2, 3}, {1, 5, 6}, {1, 6, 2}, {3, 2, 6}, {3, 6, 7}}};

// This is a special case, where since we require a convex and
// consistent ordering, we can split this into 5 tets, which is
// not the same (6 tet) discretization that could exist
// through triangulating the faces.
static constexpr std::array<std::array<UnsignedIndex_t, 4>, 5>
    tet_decomposition{
        {{0, 1, 2, 5}, {0, 2, 7, 5}, {0, 4, 5, 7}, {0, 2, 3, 7}, {2, 5, 6, 7}}};
}  // namespace hexahedron_triangulation

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType>
HexahedronSpecialization<Derived, VertexType>::generateHalfEdgeVersion(
    void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void HexahedronSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  a_half_edge_version->resize(24, 8, 6);

  for (UnsignedIndex_t v = 0; v < 8; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  static constexpr std::array<UnsignedIndex_t, 24> ending_vertex_mapping{
      {7, 6, 5, 4, 1, 2, 3, 0, 3, 7, 4, 0, 5, 6, 2, 1, 4, 5, 1, 0, 2, 6, 7, 3}};
  static constexpr std::array<UnsignedIndex_t, 24> previous_half_edge_mapping{
      {3,  0,  1,  2,  7,  4,  5,  6,  11, 8,  9,  10,
       15, 12, 13, 14, 19, 16, 17, 18, 23, 20, 21, 22}};
  static constexpr std::array<UnsignedIndex_t, 24> next_half_edge_mapping{
      {1,  2,  3,  0,  5,  6,  7,  4,  9,  10, 11, 8,
       13, 14, 15, 12, 17, 18, 19, 16, 21, 22, 23, 20}};
  static constexpr std::array<UnsignedIndex_t, 24> face_mapping{
      {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5}};
  static constexpr std::array<UnsignedIndex_t, 24> opposite_half_edge_mapping{
      {10, 22, 13, 17, 19, 15, 20, 8, 7, 23, 0, 16,
       18, 2,  21, 5,  11, 3,  12, 4, 6, 14, 1, 9}};
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

namespace hexahedron_detail {
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

}  // namespace hexahedron_detail

template <class Derived, class VertexType>
PlanarLocalizer HexahedronSpecialization<Derived, VertexType>::getLocalizer(
    void) const {
  PlanarLocalizer reconstruction_to_return;
  reconstruction_to_return.setNumberOfPlanes(6);
  Pt volume_centroid = this->calculateCentroid();

  auto hex_face = ProxyTri<Derived>::fromNoExistencePlane(
      static_cast<const Derived&>(*this), {4, 7, 5});
  reconstruction_to_return[0] =
      hexahedron_detail::planeFromTri(hex_face, volume_centroid);

  hex_face = ProxyTri<Derived>::fromNoExistencePlane(
      static_cast<const Derived&>(*this), {0, 1, 2});
  reconstruction_to_return[1] =
      hexahedron_detail::planeFromTri(hex_face, volume_centroid);

  hex_face = ProxyTri<Derived>::fromNoExistencePlane(
      static_cast<const Derived&>(*this), {4, 0, 3});
  reconstruction_to_return[2] =
      hexahedron_detail::planeFromTri(hex_face, volume_centroid);

  hex_face = ProxyTri<Derived>::fromNoExistencePlane(
      static_cast<const Derived&>(*this), {1, 5, 2});
  reconstruction_to_return[3] =
      hexahedron_detail::planeFromTri(hex_face, volume_centroid);

  hex_face = ProxyTri<Derived>::fromNoExistencePlane(
      static_cast<const Derived&>(*this), {5, 1, 0});
  reconstruction_to_return[4] =
      hexahedron_detail::planeFromTri(hex_face, volume_centroid);

  hex_face = ProxyTri<Derived>::fromNoExistencePlane(
      static_cast<const Derived&>(*this), {2, 7, 3});
  reconstruction_to_return[5] =
      hexahedron_detail::planeFromTri(hex_face, volume_centroid);

  return reconstruction_to_return;
}

template <class Derived, class VertexType>
constexpr UnsignedIndex_t HexahedronSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) {
  return static_cast<UnsignedIndex_t>(
      hexahedron_triangulation::tet_decomposition.size());
}

template <class Derived, class VertexType>
constexpr std::array<UnsignedIndex_t, 4>
HexahedronSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) {
  assert(a_tet <
         HexahedronSpecialization::getNumberOfSimplicesInDecomposition());
  return {hexahedron_triangulation::tet_decomposition[a_tet][0],
          hexahedron_triangulation::tet_decomposition[a_tet][1],
          hexahedron_triangulation::tet_decomposition[a_tet][2],
          hexahedron_triangulation::tet_decomposition[a_tet][3]};
}

template <class Derived, class VertexType>
ProxyTet<Derived>
HexahedronSpecialization<Derived, VertexType>::getSimplexFromDecomposition(
    const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYHEDRONS_HEXAHEDRON_TPP_
