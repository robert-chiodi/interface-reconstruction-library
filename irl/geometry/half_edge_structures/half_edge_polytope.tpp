// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYTOPE_TPP_
#define IRL_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYTOPE_TPP_

#include <vector>

#include "irl/data_structures/unordered_map.h"

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
  initial_half_edge_storage_size_m =
      a_number_of_half_edges * sizeof(HalfEdgeType);
  initial_vertex_storage_size_m = a_number_of_vertices * sizeof(VertexType);
  initial_face_storage_size_m = a_number_of_faces * sizeof(FaceType);
  storage_m.resize(initial_half_edge_storage_size_m +
                   initial_vertex_storage_size_m + initial_face_storage_size_m);
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
UnsignedIndex_t
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNumberOfInitialFaces(void) const {
  return static_cast<UnsignedIndex_t>(initial_face_storage_size_m /
                                      sizeof(FaceType));
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
UnsignedIndex_t
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNumberOfInitialVertices(void)
    const {
  return static_cast<UnsignedIndex_t>(initial_vertex_storage_size_m /
                                      sizeof(VertexType));
}

// INDEXING FROM ASSUMED STORAGE ORDER OF [HALF EDGES, VERTICES, FACES]
template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
HalfEdgeType &
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getHalfEdge(const UnsignedIndex_t
                                                           a_index) {
  return *reinterpret_cast<HalfEdgeType *>(
      storage_m[a_index * sizeof(HalfEdgeType)]);
}

// INDEXING FROM ASSUMED STORAGE ORDER OF [HALF EDGES, VERTICES, FACES]
template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
const HalfEdgeType &
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getHalfEdge(const UnsignedIndex_t
                                                           a_index) const {
  return *reinterpret_cast<const HalfEdgeType *>(
      storage_m[a_index * sizeof(HalfEdgeType)]);
}

// INDEXING FROM ASSUMED STORAGE ORDER OF [HALF EDGES, VERTICES, FACES]
template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
VertexType &
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getVertex(const UnsignedIndex_t
                                                         a_index) {
  return *reinterpret_cast<VertexType *>(
      storage_m[initial_half_edge_storage_size_m +
                a_index * sizeof(VertexType)]);
}

// INDEXING FROM ASSUMED STORAGE ORDER OF [HALF EDGES, VERTICES, FACES]
template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
const VertexType &
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getVertex(const UnsignedIndex_t
                                                         a_index) const {
  return *reinterpret_cast<const VertexType *>(
      storage_m[initial_half_edge_storage_size_m +
                a_index * sizeof(VertexType)]);
}

// INDEXING FROM ASSUMED STORAGE ORDER OF [HALF EDGES, VERTICES, FACES]
template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
FaceType &HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType,
                           kMaxHalfEdges, kMaxVertices,
                           kMaxFaces>::getFace(const UnsignedIndex_t a_index) {
  return *reinterpret_cast<FaceType *>(
      storage_m[initial_half_edge_storage_size_m +
                initial_vertex_storage_size_m + a_index * sizeof(FaceType)]);
}
// INDEXING FROM ASSUMED STORAGE ORDER OF [HALF EDGES, VERTICES, FACES]
template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
const FaceType &
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getFace(const UnsignedIndex_t
                                                       a_index) const {
  return *reinterpret_cast<const FaceType *>(
      storage_m[initial_half_edge_storage_size_m +
                initial_vertex_storage_size_m + a_index * sizeof(FaceType)]);
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
HalfEdgeType *
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNewHalfEdge(void) {
  auto new_object = storage_m.template getNewObject<HalfEdgeType>();
  *new_object = HalfEdgeType();
  return new_object;
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
HalfEdgeType *
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNewHalfEdge(const HalfEdgeType
                                                              &a_half_edge) {
  auto new_object = storage_m.template getNewObject<HalfEdgeType>();
  *new_object = a_half_edge;
  return new_object;
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
HalfEdgeType *
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNewHalfEdge(HalfEdgeType
                                                              &&a_half_edge) {
  auto new_object = storage_m.template getNewObject<HalfEdgeType>();
  *new_object = std::move(a_half_edge);
  return new_object;
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
VertexType *
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNewVertex(void) {
  auto new_object = storage_m.template getNewObject<VertexType>();
  *new_object = VertexType();
  return new_object;
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
VertexType *
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNewVertex(VertexType &&a_vertex) {

  auto new_object = storage_m.template getNewObject<VertexType>();
  *new_object = std::move(a_vertex);
  return new_object;
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
FaceType *
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNewFace(void) {
  auto new_object = storage_m.template getNewObject<FaceType>();
  *new_object = FaceType();
  return new_object;
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
FaceType *
HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                 kMaxVertices, kMaxFaces>::getNewFace(FaceType &&a_face) {
  auto new_object = storage_m.template getNewObject<FaceType>();
  *new_object = std::move(a_face);
  return new_object;
}

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
template <class GeometryType>
void HalfEdgePolytope<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::setVertexLocations(const GeometryType &a_geometry) {
  assert(a_geometry.getNumberOfVertices() ==
         this->getNumberOfInitialVertices());
  auto vertex = reinterpret_cast<VertexType *>(
      storage_m[initial_half_edge_storage_size_m]);
  for (UnsignedIndex_t v = 0; v < a_geometry.getNumberOfVertices(); ++v) {
    (vertex++)->setLocation(a_geometry[v]);
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

// template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
//           UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
//           UnsignedIndex_t kMaxFaces>
// inline std::ostream &
// operator<<(std::ostream &out,
//            const HalfEdgePolytope<PtType, VertexType, HalfEdgeType, FaceType,
//                                   kMaxHalfEdges, kMaxVertices, kMaxFaces>
//                &a_polytope) {
//   // Mapping of VertexPointer to unique id.
//   unordered_map<const VertexType *, UnsignedIndex_t> unique_vertices;
//   for (UnsignedIndex_t n = 0; n < a_polytope.getNumberOfVertices(); ++n) {
//     unique_vertices[&a_polytope.getVertex(n)] = n;
//   }
//   assert(unique_vertices.size() == a_polytope.getNumberOfVertices());

//   // VTK XML Header
//   out << "<?xml version=\"1.0\"?>\n";
//   out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
//          "byte_order=\"LittleEndian\">\n";
//   out << "<UnstructuredGrid>\n";
//   out << "<Piece NumberOfPoints=\"" << a_polytope.getNumberOfVertices()
//       << "\" NumberOfCells=\"" << a_polytope.getNumberOfFaces() << "\">\n";

//   // Vertex Locations
//   out << "<Points>\n";
//   out << "<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
//   for (UnsignedIndex_t n = 0; n < a_polytope.getNumberOfVertices(); ++n) {
//     const auto &vert_pt = a_polytope.getVertex(n).getLocation();
//     out << vert_pt[0] << " " << vert_pt[1] << " " << vert_pt[2] << '\n';
//   }
//   out << "</DataArray>\n";
//   out << "</Points>\n";

//   // Creation of "Cell" from group of polygons
//   out << "<Cells>\n";
//   // Connectivity of vertices
//   out << "<DataArray type=\"Int32\" Name=\"connectivity\" "
//          "format=\"ascii\">\n";
//   std::vector<int> face_sizes(a_polytope.getNumberOfFaces());
//   for (UnsignedIndex_t n = 0; n < a_polytope.getNumberOfFaces(); ++n) {
//     const auto &face = a_polytope.getFace(n);
//     int current_face_size = 0;
//     auto current_half_edge = face.getStartingHalfEdge();
//     do {
//       ++current_face_size;
//       out << unique_vertices[current_half_edge->getVertex()] << " ";
//       current_half_edge = current_half_edge->getNextHalfEdge();
//     } while (current_half_edge != face.getStartingHalfEdge());
//     face_sizes[n] = current_face_size;
//     out << '\n';
//   }
//   out << "</DataArray>\n";

//   // Starting offset for connectivity of each polygon in the cell
//   out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
//   int offset = 0;
//   for (const auto &element : face_sizes) {
//     offset += element;
//     out << offset << " ";
//   }
//   out << '\n';
//   out << "</DataArray>\n";

//   // Cell type - General Polygon type
//   out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
//   for (UnsignedIndex_t n = 0; n < a_polytope.getNumberOfFaces(); ++n) {
//     out << "7 \n"; // General polygon type
//   }
//   out << "</DataArray>\n";
//   out << "</Cells>\n";
//   out << "</Piece>\n";
//   out << "</UnstructuredGrid>\n";
//   out << "</VTKFile>\n";

//   return out;
// }

} // namespace IRL

#endif // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_HALF_EDGE_POLYTOPE_TPP_
