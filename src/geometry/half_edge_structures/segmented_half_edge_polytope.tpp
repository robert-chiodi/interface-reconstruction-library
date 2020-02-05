// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYTOPE_TPP_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYTOPE_TPP_

#include "src/helpers/mymath.h"

namespace IRL {

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
void SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces, kMaxVertices>::
    setNumberOfFaces(const UnsignedIndex_t a_number_of_faces) {
  faces_m.resize(a_number_of_faces);
  this->checkIfStaticAllocationExceeded();
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
void SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces, kMaxVertices>::
    setNumberOfVertices(const UnsignedIndex_t a_size) {
  vertices_m.resize(a_size);
  this->checkIfStaticAllocationExceeded();
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
UnsignedIndex_t SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                                          kMaxVertices>::getNumberOfFaces(void)
    const {
  return static_cast<UnsignedIndex_t>(faces_m.size());
}
template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
UnsignedIndex_t
SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                          kMaxVertices>::getNumberOfVertices(void) const {
  return static_cast<UnsignedIndex_t>(vertices_m.size());
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
FaceType *SegmentedHalfEdgePolytope<
    FaceType, VertexType, kMaxFaces,
    kMaxVertices>::operator[](const UnsignedIndex_t a_index) {
  return faces_m[a_index];
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
const FaceType *SegmentedHalfEdgePolytope<
    FaceType, VertexType, kMaxFaces,
    kMaxVertices>::operator[](const UnsignedIndex_t a_index) const {
  return faces_m[a_index];
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
void SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                               kMaxVertices>::addFace(FaceType *a_face) {
  faces_m.push_back(a_face);
  this->checkIfStaticAllocationExceeded();
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
FaceType *&SegmentedHalfEdgePolytope<
    FaceType, VertexType, kMaxFaces,
    kMaxVertices>::getFacePointer(const UnsignedIndex_t a_index) {
  assert(a_index < faces_m.size());
  return faces_m[a_index];
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
void SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                               kMaxVertices>::addVertex(VertexType *a_vertex) {
  vertices_m.push_back(a_vertex);
  this->checkIfStaticAllocationExceeded();
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
VertexType *SegmentedHalfEdgePolytope<
    FaceType, VertexType, kMaxFaces,
    kMaxVertices>::getVertex(const UnsignedIndex_t a_index) {
  assert(a_index < static_cast<UnsignedIndex_t>(vertices_m.size()));
  return vertices_m[a_index];
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
const VertexType *SegmentedHalfEdgePolytope<
    FaceType, VertexType, kMaxFaces,
    kMaxVertices>::getVertex(const UnsignedIndex_t a_index) const {
  assert(a_index < static_cast<UnsignedIndex_t>(vertices_m.size()));
  return vertices_m[a_index];
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
VertexType *&SegmentedHalfEdgePolytope<
    FaceType, VertexType, kMaxFaces,
    kMaxVertices>::getVertexPointer(const UnsignedIndex_t a_index) {
  assert(a_index < vertices_m.size());
  return vertices_m[a_index];
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
Pt SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                             kMaxVertices>::getLowerLimits(void) const {
  Pt pt_to_return(DBL_MAX, DBL_MAX, DBL_MAX);
  for (const auto &vertex : vertices_m) {
    const auto &pt = vertex->getLocation();
    pt_to_return[0] = std::min(pt_to_return[0], pt[0]);
    pt_to_return[1] = std::min(pt_to_return[1], pt[1]);
    pt_to_return[2] = std::min(pt_to_return[2], pt[2]);
  }
  return pt_to_return;
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
Pt SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                             kMaxVertices>::getUpperLimits(void) const {
  Pt pt_to_return(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for (const auto &vertex : vertices_m) {
    const auto &pt = vertex->getLocation();
    pt_to_return[0] = std::max(pt_to_return[0], pt[0]);
    pt_to_return[1] = std::max(pt_to_return[1], pt[1]);
    pt_to_return[2] = std::max(pt_to_return[2], pt[2]);
  }
  return pt_to_return;
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
std::array<Pt, 2> SegmentedHalfEdgePolytope<
    FaceType, VertexType, kMaxFaces, kMaxVertices>::getBoundingBox(void) const {
  std::array<Pt, 2> bounding_box{
      {Pt(DBL_MAX, DBL_MAX, DBL_MAX), Pt(-DBL_MAX, -DBL_MAX, -DBL_MAX)}};
  for (const auto &vertex : vertices_m) {
    const auto &pt = vertex->getLocation();
    for (UnsignedIndex_t dim = 0; dim < 3; ++dim) {
      bounding_box[0][dim] = std::min(bounding_box[0][dim], pt[dim]);
      bounding_box[1][dim] = std::max(bounding_box[1][dim], pt[dim]);
    }
  }
  return bounding_box;
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
int SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces, kMaxVertices>::
    calculateAndStoreDistanceToVertices(const Plane &a_plane) {
  // No actual object
  if (this->getNumberOfVertices() == 0) {
    return -1;
  }

  // Plane is set to non-intersect
  if (squaredMagnitude(a_plane.normal()) < DBL_MIN) {
    return a_plane.distance() > 0.0 ? -1 : 1;
  }

  // Leaving this here in case it later becomes useful.
  // Calculates a bounding box and first checks for intersections on it.
  // If intersection exists, go on to check actual object.
  // Currently, this actually is more expensive than just always checking
  // for an intersection using the vertices.
  //  static constexpr UnsignedIndex_t small_number_of_vertices = 8;
  //  if(this->getNumberOfVertices() > small_number_of_vertices){
  //    auto bounding_box = this->getBoundingBox();
  //    for(UnsignedIndex_t dim = 0; dim < 3; ++dim){
  //      if(a_plane.normal()[dim] < 0.0){
  //        auto tmp = bounding_box[0][dim];
  //      bounding_box[0][dim] = bounding_box[1][dim];
  //      bounding_box[1][dim] = tmp;
  //    }
  //    }
  //    std::array<double, 2>
  //    bounding_box_distances{{a_plane.signedDistanceToPoint(bounding_box[0]),
  //      a_plane.signedDistanceToPoint(bounding_box[1])}};
  //    if(bounding_box_distances[0]*bounding_box_distances[1] >= 0.0){
  //      if(bounding_box_distances[0] != 0.0){
  //        return bounding_box_distances[0] < 0.0 ? -1 : 1;
  //      } else {
  //        return bounding_box_distances[1] < 0.0 ? -1 : 1;
  //      }
  //    }
  //  }

  const auto nvert = vertices_m.size();
  for (UnsignedIndex_t v = 0; v < nvert; ++v) {
    auto &vertex = *vertices_m[v];
    vertex.calculateDistanceToPlane(a_plane);
    vertex.setClip(vertex.getDistance() > 0.0);
  }

  // -1 is fully under
  // 1 is fully above
  // 0 is intersected
  if (vertices_m[0]->isClipped()) {
    for (UnsignedIndex_t v = 1; v < nvert; ++v) {
      const auto &vertex = *vertices_m[v];
      if (vertex.isNotClipped()) {
        return 0;
      }
    }
    return 1;
  } else {
    for (UnsignedIndex_t v = 1; v < nvert; ++v) {
      const auto &vertex = *vertices_m[v];
      if (vertex.isClipped()) {
        return 0;
      }
    }
    return -1;
  }
  return 0;  // Should never reach
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
template <class HalfEdgePolytopeType>
void SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces, kMaxVertices>::
    clear(HalfEdgePolytopeType *a_complete_polytope) {
  for (auto &vertex : vertices_m) {
    a_complete_polytope->freeVertexAndHalfEdges(vertex);
  }
  for (auto &face : faces_m) {
    a_complete_polytope->freeFace(face);
  }
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
void SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces, kMaxVertices>::
    removeFace(const UnsignedIndex_t a_index) {
  faces_m.erase(faces_m.begin() + a_index);
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
void SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces, kMaxVertices>::
    removeVertex(const UnsignedIndex_t a_index) {
  vertices_m.erase(vertices_m.begin() + a_index);
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
typename SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                                   kMaxVertices>::iterator
SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces, kMaxVertices>::begin(
    void) noexcept {
  return faces_m.begin();
}
template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
typename SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                                   kMaxVertices>::const_iterator
SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces, kMaxVertices>::begin(
    void) const noexcept {
  return this->cbegin();
}
template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
typename SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                                   kMaxVertices>::const_iterator
SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                          kMaxVertices>::cbegin(void) const noexcept {
  return faces_m.cbegin();
}
template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
typename SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                                   kMaxVertices>::iterator
SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces, kMaxVertices>::end(
    void) noexcept {
  return faces_m.end();
}
template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
typename SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                                   kMaxVertices>::const_iterator
SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces, kMaxVertices>::end(
    void) const noexcept {
  return this->cend();
}
template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
typename SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                                   kMaxVertices>::const_iterator
SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces, kMaxVertices>::cend(
    void) const noexcept {
  return faces_m.cend();
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
void SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces, kMaxVertices>::
    checkIfStaticAllocationExceeded(void) const {
#ifndef NDEBUG_PERF
  if (faces_m.capacity() > kMaxFaces) {
    std::cout << "Static allocation size for SmallVector exceeded in "
                 "SegmentedHalfEdgePolytope for faces_m. Expect performance "
                 "penalty if this happens frequently."
              << std::endl;
  }
  if (vertices_m.capacity() > kMaxVertices) {
    std::cout << "Static allocation size for SmallVector exceeded in "
                 "SegmentedHalfEdgePolytope for vertices_m. Expect "
                 "performance penalty if this happens frequently."
              << std::endl;
  }
#endif
}

template <class FaceType, class VertexType, UnsignedIndex_t kMaxFaces,
          UnsignedIndex_t kMaxVertices>
inline std::ostream &operator<<(
    std::ostream &out,
    const SegmentedHalfEdgePolytope<FaceType, VertexType, kMaxFaces,
                                    kMaxVertices> &a_polytope) {
  UnsignedIndex_t global_number_of_verts = 0;
  for (auto &face : a_polytope) {
    if (face == &getOpenBoundaryFace<FaceType>()) {
      continue;
    }
    auto current_half_edge = face->getStartingHalfEdge();
    do {
      ++global_number_of_verts;
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != face->getStartingHalfEdge());
  }

  out << "<?xml version=\"1.0\"?>" << std::endl;
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
         "byte_order=\"LittleEndian\">"
      << std::endl;
  out << "<UnstructuredGrid>" << std::endl;
  out << "<Piece NumberOfPoints=\"" << global_number_of_verts
      << "\" NumberOfCells=\"" << a_polytope.getNumberOfFaces() << "\">"
      << std::endl;
  //////// Vertices
  out << "<Points>" << std::endl;
  out << "<DataArray type=\"Float32\" NumberOfComponents=\"3\">" << std::endl;
  for (auto &face : a_polytope) {
    if (face == &getOpenBoundaryFace<FaceType>()) {
      continue;
    }
    auto current_half_edge = face->getStartingHalfEdge();
    do {
      const auto &vert_pt = current_half_edge->getVertex()->getLocation();
      out << vert_pt[0] << " " << vert_pt[1] << " " << vert_pt[2] << std::endl;
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != face->getStartingHalfEdge());
  }
  out << "</DataArray>" << std::endl;
  out << "</Points>" << std::endl;
  //////// Cells
  out << "<Cells>" << std::endl;
  // Connectivity
  out << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
      << std::endl;
  UnsignedIndex_t number_of_verts = 0;
  for (auto &face : a_polytope) {
    if (face == &getOpenBoundaryFace<FaceType>()) {
      continue;
    }
    auto current_half_edge = face->getStartingHalfEdge();
    do {
      out << number_of_verts << " ";
      ++number_of_verts;
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != face->getStartingHalfEdge());
    out << std::endl;
  }
  out << "</DataArray>" << std::endl;
  // Offsets
  out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
      << std::endl;
  UnsignedIndex_t offset = 0;
  for (auto &face : a_polytope) {
    if (face == &getOpenBoundaryFace<FaceType>()) {
      continue;
    }
    auto current_half_edge = face->getStartingHalfEdge();
    do {
      ++offset;
      current_half_edge = current_half_edge->getNextHalfEdge();
    } while (current_half_edge != face->getStartingHalfEdge());
    out << offset << " ";
  }
  out << std::endl;
  out << "</DataArray>" << std::endl;
  // Cell type
  out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
      << std::endl;
  for (auto &face : a_polytope) {
    if (face == &getOpenBoundaryFace<FaceType>()) {
      continue;
    }
    out << "7" << std::endl;  // General polygon type
  }
  out << "</DataArray>" << std::endl;
  out << "</Cells>" << std::endl;
  out << "</Piece>" << std::endl;
  out << "</UnstructuredGrid>" << std::endl;
  out << "</VTKFile>" << std::endl;

  return out;
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_SEGMENTED_HALF_EDGE_POLYTOPE_TPP_