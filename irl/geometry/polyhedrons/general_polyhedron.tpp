// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_POLYHEDRONS_GENERAL_POLYHEDRON_TPP_
#define IRL_GEOMETRY_POLYHEDRONS_GENERAL_POLYHEDRON_TPP_

#include <cassert>

#include "irl/geometry/half_edge_structures/half_edge.h"

namespace IRL {

template <class Derived, class VertexType>
GeneralPolyhedronSpecialization<
    Derived, VertexType>::GeneralPolyhedronSpecialization(void)
    : connectivity_m(nullptr) {}

template <class Derived, class VertexType>
GeneralPolyhedronSpecialization<Derived, VertexType>::
    GeneralPolyhedronSpecialization(
        const PolyhedronConnectivity* a_connectivity)
    : connectivity_m{a_connectivity} {}

template <class Derived, class VertexType>
HalfEdgePolyhedron<VertexType>
GeneralPolyhedronSpecialization<Derived, VertexType>::generateHalfEdgeVersion(
    void) const {
  HalfEdgePolyhedron<VertexType> half_edge_version;
  this->setHalfEdgeVersion(&half_edge_version);
  return half_edge_version;
}

template <class Derived, class VertexType>
template <class HalfEdgePolyhedronType>
void GeneralPolyhedronSpecialization<Derived, VertexType>::setHalfEdgeVersion(
    HalfEdgePolyhedronType* a_half_edge_version) const {
  using HalfEdgeType = typename HalfEdgePolyhedronType::half_edge_type;

  assert(connectivity_m != nullptr);

  const auto number_of_half_edges = static_cast<UnsignedIndex_t>(
      connectivity_m->ending_vertex_mapping.size());
  const auto number_of_vertices = connectivity_m->number_of_vertices;
  assert(number_of_half_edges >= 12);
  assert(number_of_vertices > 3);
  const auto number_of_faces =
      number_of_half_edges / 2 - number_of_vertices + 2;
  a_half_edge_version->resize(number_of_half_edges, number_of_vertices,
                              number_of_faces);

  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    a_half_edge_version->getVertex(v).setLocation((*this)[v]);
  }

  for (UnsignedIndex_t n = 0; n < number_of_half_edges; ++n) {
    HalfEdgeType& current_half_edge = a_half_edge_version->getHalfEdge(n);
    current_half_edge = HalfEdgeType(
        &a_half_edge_version->getVertex(
            connectivity_m->ending_vertex_mapping.at(n)),
        &a_half_edge_version->getHalfEdge(
            connectivity_m->previous_half_edge_mapping.at(n)),
        &a_half_edge_version->getHalfEdge(
            connectivity_m->next_half_edge_mapping.at(n)),
        &a_half_edge_version->getFace(connectivity_m->face_mapping.at(n)));
    current_half_edge.setOppositeHalfEdge(&a_half_edge_version->getHalfEdge(
        connectivity_m->opposite_half_edge_mapping.at(n)));
    current_half_edge.getFace()->setStartingHalfEdge(&current_half_edge);
    current_half_edge.getVertex()->setHalfEdge(&current_half_edge);
  }
}

template <class Derived, class VertexType>
std::array<UnsignedIndex_t, 4>
GeneralPolyhedronSpecialization<Derived, VertexType>::
    getSimplexIndicesFromDecomposition(const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  assert(connectivity_m != nullptr);
  return {connectivity_m->face_triangle_decomposition[a_tet][0],
          connectivity_m->face_triangle_decomposition[a_tet][1],
          connectivity_m->face_triangle_decomposition[a_tet][2],
          connectivity_m->datum_index};
}

template <class Derived, class VertexType>
UnsignedIndex_t GeneralPolyhedronSpecialization<
    Derived, VertexType>::getNumberOfSimplicesInDecomposition(void) const {
  assert(connectivity_m != nullptr);
  return static_cast<UnsignedIndex_t>(
      connectivity_m->face_triangle_decomposition.size());
}

template <class Derived, class VertexType>
ProxyTet<Derived> GeneralPolyhedronSpecialization<Derived, VertexType>::
    getSimplexFromDecomposition(const UnsignedIndex_t a_tet) const {
  assert(a_tet < this->getNumberOfSimplicesInDecomposition());
  return {static_cast<const Derived&>(*this),
          this->getSimplexIndicesFromDecomposition(a_tet)};
}

template <class Derived, class VertexType>
void GeneralPolyhedronSpecialization<Derived, VertexType>::resetConnectivity(
    const PolyhedronConnectivity* a_connectivity) {
  connectivity_m = a_connectivity;
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
template <class VertexListType>
StoredGeneralPolyhedron<VertexType, kStaticAllocSize>::StoredGeneralPolyhedron(
    const VertexListType& a_vertex_list,
    const PolyhedronConnectivity* a_connectivity)
    : StoredGeneralPolyhedron<VertexType, kStaticAllocSize>::Base(
          a_connectivity) {
  assert(a_vertex_list.size() == this->connectivity_m->number_of_vertices);
  vertex_list_m.setNumberOfPts(static_cast<UnsignedIndex_t>(a_vertex_list.size()));
  std::copy(a_vertex_list.begin(), a_vertex_list.end(), vertex_list_m.begin());
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
StoredGeneralPolyhedron<VertexType, kStaticAllocSize>
StoredGeneralPolyhedron<VertexType, kStaticAllocSize>::fromRawPtPointer(
    const UnsignedIndex_t a_number_of_pts, const VertexType* a_array_of_pts,
    const PolyhedronConnectivity* a_connectivity) {
  return StoredGeneralPolyhedron(a_number_of_pts, a_array_of_pts,
                                 a_connectivity);
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
StoredGeneralPolyhedron<VertexType, kStaticAllocSize>
StoredGeneralPolyhedron<VertexType, kStaticAllocSize>::fromRawDoublePointer(
    const UnsignedIndex_t a_number_of_pts, const double* a_array_of_locs,
    const PolyhedronConnectivity* a_connectivity) {
  return StoredGeneralPolyhedron(a_number_of_pts, a_array_of_locs,
                                 a_connectivity);
}

template <class VertexType, UnsignedIndex_t kStaticAllocSize>
VertexType& StoredGeneralPolyhedron<VertexType, kStaticAllocSize>::access(
    const UnsignedIndex_t a_index) {
  assert(a_index < this->getNumberOfVerticesInObject());
  return vertex_list_m[a_index];
}
template <class VertexType, UnsignedIndex_t kStaticAllocSize>
const VertexType& StoredGeneralPolyhedron<VertexType, kStaticAllocSize>::access(
    const UnsignedIndex_t a_index) const {
  assert(a_index < this->getNumberOfVerticesInObject());
  return vertex_list_m[a_index];
}
template <class VertexType, UnsignedIndex_t kStaticAllocSize>
UnsignedIndex_t StoredGeneralPolyhedron<
    VertexType, kStaticAllocSize>::getNumberOfVerticesInObject(void) const {
  return static_cast<UnsignedIndex_t>(vertex_list_m.getNumberOfPts());
}
template <class VertexType, UnsignedIndex_t kStaticAllocSize>
void StoredGeneralPolyhedron<VertexType, kStaticAllocSize>::setNumberOfVertices(
    const UnsignedIndex_t a_number) {
  vertex_list_m.setNumberOfPts(a_number);
}
template <class VertexType, UnsignedIndex_t kStaticAllocSize>
StoredGeneralPolyhedron<VertexType, kStaticAllocSize>::StoredGeneralPolyhedron(
    const UnsignedIndex_t a_number_of_pts, const VertexType* a_array_of_pts,
    const PolyhedronConnectivity* a_connectivity)
    : StoredGeneralPolyhedron<VertexType, kStaticAllocSize>::Base(
          a_connectivity) {
  assert(a_number_of_pts == this->connectivity_m->number_of_vertices);
  vertex_list_m.setNumberOfPts(a_number_of_pts);
  std::copy(a_array_of_pts, a_array_of_pts + a_number_of_pts,
            vertex_list_m.size());
}
template <class VertexType, UnsignedIndex_t kStaticAllocSize>
StoredGeneralPolyhedron<VertexType, kStaticAllocSize>::StoredGeneralPolyhedron(
    const UnsignedIndex_t a_number_of_pts, const double* a_array_of_locs,
    const PolyhedronConnectivity* a_connectivity)
    : StoredGeneralPolyhedron<VertexType, kStaticAllocSize>::Base(
          a_connectivity) {
  assert(a_number_of_pts == this->connectivity_m->number_of_vertices);
  vertex_list_m.setNumberOfPts(a_number_of_pts);
  for (UnsignedIndex_t n = 0; n < a_number_of_pts; ++n) {
    for (UnsignedIndex_t d = 0; d < a_number_of_pts; ++d) {
      vertex_list_m[n][d] = a_array_of_locs[n * 3 + d];
    }
  }
}

}  // namespace IRL

#endif // IRL_GEOMETRY_POLYHEDRONS_GENERAL_POLYHEDRON_TPP_
