// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_POLYHEDRONS_POLYHEDRON_CONNECTIVITY_H_
#define IRL_GEOMETRY_POLYHEDRONS_POLYHEDRON_CONNECTIVITY_H_

#include <array>
#include <initializer_list>

#include "irl/parameters/defined_types.h"

namespace IRL {

class PolyhedronConnectivity {
 public:
  PolyhedronConnectivity(void) = default;

  // FaceBREPType should be a 2D array type that enumerates
  // the vertex numbers for each face. Both array types must
  // return their size with a `size()` method, and allow access
  // via `operator[]`. Examples include
  // std::vector<std::vector<UnsignedIndex_t>>,
  // std::vector<std::array<UnsignedIndex_t, number>>, or other such
  // combinations.
  template <class FaceBREPType>
  PolyhedronConnectivity(const FaceBREPType& a_face_brep);

  ~PolyhedronConnectivity(void) = default;

 public:
  UnsignedIndex_t number_of_vertices;
  UnsignedIndex_t datum_index;
  std::vector<std::array<UnsignedIndex_t, 3>> face_triangle_decomposition;
  std::vector<UnsignedIndex_t> ending_vertex_mapping;
  std::vector<UnsignedIndex_t> previous_half_edge_mapping;
  std::vector<UnsignedIndex_t> next_half_edge_mapping;
  std::vector<UnsignedIndex_t> face_mapping;
  std::vector<UnsignedIndex_t> opposite_half_edge_mapping;
};

}  // namespace IRL

#include "irl/geometry/polyhedrons/polyhedron_connectivity.tpp"
#endif // IRL_GEOMETRY_POLYHEDRONS_POLYHEDRON_CONNECTIVITY_H_
