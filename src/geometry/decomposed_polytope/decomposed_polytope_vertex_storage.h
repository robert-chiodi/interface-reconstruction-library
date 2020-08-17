// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_DECOMPOSED_POLYTOPE_DECOMPOSED_POLYTOPE_VERTEX_STORAGE_H_
#define SRC_GEOMETRY_DECOMPOSED_POLYTOPE_DECOMPOSED_POLYTOPE_VERTEX_STORAGE_H_

#include <vector>

#include "src/data_structures/chained_block_storage.h"
#include "src/data_structures/small_vector.h"

namespace IRL {

template <class VertexType>
class DecomposedPolytopeVertexStorage {
 public:
  using pt_type = VertexType;

  template <class GeometryType>
  void resetFromGeometry(const GeometryType& a_geometry) {
    vertex_storage_m.resize(a_geometry.getNumberOfVertices());
    distance_m.resize(a_geometry.getNumberOfVertices());
    for (UnsignedIndex_t n = 0; n < a_geometry.getNumberOfVertices(); ++n) {
      vertex_storage_m[n] = pt_type(a_geometry[n]);
    }
  }

  pt_type& operator[](const UnsignedIndex_t a_index) {
    return vertex_storage_m[a_index];
  }

  const pt_type& operator[](const UnsignedIndex_t a_index) const {
    return vertex_storage_m[a_index];
  }

  void push_back(const VertexType a_vertex) {
    vertex_storage_m.push_back(a_vertex);
    distance_m.push_back(0.0);
  }

  UnsignedIndex_t size(void) const {
    return static_cast<UnsignedIndex_t>(vertex_storage_m.size());
  }
  void resize(const UnsignedIndex_t a_new_size) {
    vertex_storage_m.resize(a_new_size);
    distance_m.resize(a_new_size);
  }

  double getDistance(const UnsignedIndex_t a_index) const {
    return distance_m[a_index];
  }

  double setDistance(const UnsignedIndex_t a_index, const double a_distance) {
    return distance_m[a_index] = a_distance;
  }

 private:
  std::vector<VertexType> vertex_storage_m;
  std::vector<double> distance_m;
};

template <class VertexType>
class DecomposedPolyhedronVertexStorage
    : public DecomposedPolytopeVertexStorage<VertexType> {
 public:
 private:
};

template <class VertexType>
class DecomposedPolygonVertexStorage
    : public DecomposedPolytopeVertexStorage<VertexType> {
 public:
  template <class GeometryType>
  void resetFromGeometry(const GeometryType& a_geometry) {
    DecomposedPolytopeVertexStorage<VertexType>::resetFromGeometry(a_geometry);
  }

  template <class BasePtType>
  void resetFromGeometry(
      const ExpandableDividedPolygon<BasePtType>& a_geometry) {
    // Strip the mask that skips the first index, allowing copying of the
    // centroid which is needed for creation of simplices.
    DecomposedPolytopeVertexStorage<VertexType>::resetFromGeometry(
        static_cast<const MaskStripper<ExpandableDividedPolygon<BasePtType>>&>(
            a_geometry));
    // Add last point that was not added due to number of vertices being 1 less
    // than correct once stripped.
    this->push_back(
        VertexType(a_geometry[a_geometry.getNumberOfVertices() - 1]));
  }

  void setPlaneOfExistence(const Plane& a_plane) {
    plane_of_existence_m = a_plane;
  }
  const Plane& getPlaneOfExistence(void) const { return plane_of_existence_m; }

 private:
  Plane plane_of_existence_m;
};

}  // namespace IRL
#endif  // SRC_GEOMETRY_DECOMPOSED_POLYTOPE_DECOMPOSED_POLYTOPE_VERTEX_STORAGE_H_
