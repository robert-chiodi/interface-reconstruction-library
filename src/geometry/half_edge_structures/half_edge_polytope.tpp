// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYTOPE_TPP_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYTOPE_TPP_

namespace IRL {

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>
HalfEdgePolytope<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::fromKnownSizes(const UnsignedIndex_t a_number_of_half_edges,
                               const UnsignedIndex_t a_number_of_vertices,
                               const UnsignedIndex_t a_number_of_faces) {
  return HalfEdgePolytope(a_number_of_half_edges, a_number_of_vertices,
                          a_number_of_faces);
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
void HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                      kMaxVertices, kMaxFaces>::reset(void) {
  this->resize(0, 0, 0);
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
void HalfEdgePolytope<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::resize(const UnsignedIndex_t a_number_of_half_edges,
                       const UnsignedIndex_t a_number_of_vertices,
                       const UnsignedIndex_t a_number_of_faces) {
  half_edges_m.resize(a_number_of_half_edges);
  vertices_m.resize(a_number_of_vertices);
  faces_m.resize(a_number_of_faces);
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
UnsignedIndex_t
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNumberOfFaces(void) const {
  return faces_m.size();
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
UnsignedIndex_t
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNumberOfVertices(void) const {
  return vertices_m.size();
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
HalfEdgeType& HalfEdgePolytope<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::getHalfEdge(const UnsignedIndex_t a_index) {
  return half_edges_m[a_index];
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
const HalfEdgeType& HalfEdgePolytope<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::getHalfEdge(const UnsignedIndex_t a_index) const {
  return half_edges_m[a_index];
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
VertexType& HalfEdgePolytope<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::getVertex(const UnsignedIndex_t a_index) {
  return vertices_m[a_index];
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
const VertexType& HalfEdgePolytope<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::getVertex(const UnsignedIndex_t a_index) const {
  return vertices_m[a_index];
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
FaceType& HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType,
                           kMaxHalfEdges, kMaxVertices,
                           kMaxFaces>::getFace(const UnsignedIndex_t a_index) {
  return faces_m[a_index];
}
template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
const FaceType& HalfEdgePolytope<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::getFace(const UnsignedIndex_t a_index) const {
  return faces_m[a_index];
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
HalfEdgeType*
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNewHalfEdge(void) {
  return &half_edges_m.getNextElement();
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
HalfEdgeType* HalfEdgePolytope<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::getNewHalfEdge(const HalfEdgeType& a_half_edge) {
  return &half_edges_m.getNextElement(a_half_edge);
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
HalfEdgeType* HalfEdgePolytope<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::getNewHalfEdge(HalfEdgeType&& a_half_edge) {
  return &half_edges_m.getNextElement(std::move(a_half_edge));
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
VertexType*
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNewVertex(void) {
  return &vertices_m.getNextElement();
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
VertexType*
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNewVertex(VertexType&& a_vertex) {
  return &vertices_m.getNextElement(std::move(a_vertex));
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
FaceType*
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNewFace(void) {
  return &faces_m.getNextElement();
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
FaceType*
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNewFace(FaceType&& a_face) {
  return &faces_m.getNextElement(std::move(a_face));
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
template <class GeometryType>
void HalfEdgePolytope<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::setVertexLocations(const GeometryType& a_geometry) {
  assert(a_geometry.getNumberOfVertices() == this->getNumberOfVertices());
  for (UnsignedIndex_t v = 0; v < a_geometry.getNumberOfVertices(); ++v) {
    vertices_m[v].setLocation(a_geometry[v]);
  }
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
HalfEdgePolytope<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::HalfEdgePolytope(const UnsignedIndex_t a_number_of_half_edges,
                                 const UnsignedIndex_t a_number_of_vertices,
                                 const UnsignedIndex_t a_number_of_faces) {
  this->resize(a_number_of_half_edges, a_number_of_vertices, a_number_of_faces);
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
inline std::ostream& operator<<(
    std::ostream& out,
    const HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType,
                           kMaxHalfEdges, kMaxVertices, kMaxFaces>&
        a_polyhedron) {
  out << "nfaces = " << a_polyhedron.getNumberOfFaces() << ";\n";
  out << "nverts = " << a_polyhedron.getNumberOfVertices() << ";\n";

  UnsignedIndex_t current_face = 1;
  for (UnsignedIndex_t f = 0; f < a_polyhedron.getNumberOfFaces(); ++f) {
    const auto& face = a_polyhedron.getFace(f);
    UnsignedIndex_t number_of_verts = 0;
    auto current_half_edge = face.getStartingHalfEdge();
    do {
      ++number_of_verts;
      const auto& vert_pt = current_half_edge->getVertex()->getLocation();
      out << "vert_on_face(1:3, " << number_of_verts << ", " << current_face
          << ") = ";
      out << "[ " << vert_pt[0] << ", " << vert_pt[1] << ", " << vert_pt[2]
          << "];\n";
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != face.getStartingHalfEdge());
    out << "nvert_for_face(" << current_face << ") = " << number_of_verts
        << "; \n";
    ++current_face;
  }
  return out;
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYTOPE_TPP_
