// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_POLYHEDRON_CONNECTIVITY_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_POLYHEDRON_CONNECTIVITY_TPP_

#include <algorithm>
#include <cassert>

namespace IRL {

template <class FaceBREPType>
PolyhedronConnectivity::PolyhedronConnectivity(
    const FaceBREPType& a_face_brep) {
  number_of_vertices = 0;
  UnsignedIndex_t number_of_half_edges = 0;
  for (const auto& face : a_face_brep) {
    for (const auto& vertex : face) {
      number_of_vertices = std::max(number_of_vertices, vertex);
    }
    number_of_half_edges += static_cast<UnsignedIndex_t>(face.size());
  }
  ++number_of_vertices;  // Since 0-indexed, size is one greater
  assert(number_of_half_edges / 2 ==
         a_face_brep.size() + number_of_vertices - 2);

  ending_vertex_mapping.resize(number_of_half_edges);
  previous_half_edge_mapping.resize(number_of_half_edges);
  next_half_edge_mapping.resize(number_of_half_edges);
  face_mapping.resize(number_of_half_edges);
  opposite_half_edge_mapping.resize(number_of_half_edges);

  std::vector<std::vector<UnsignedIndex_t>> half_edge_mapping(
      number_of_vertices, std::vector<UnsignedIndex_t>(number_of_vertices));

  std::vector<UnsignedIndex_t> vertex_visited(number_of_vertices, 0);

  // Generate connectivity for all half edges except opposite
  UnsignedIndex_t next_half_edge = 0;
  for (UnsignedIndex_t f = 0; f < a_face_brep.size(); ++f) {
    const auto& face = a_face_brep[f];

    ending_vertex_mapping[next_half_edge] = face[0];
    previous_half_edge_mapping[next_half_edge] =
        next_half_edge + static_cast<UnsignedIndex_t>(face.size()) - 1;
    next_half_edge_mapping[next_half_edge] = next_half_edge + 1;
    face_mapping[next_half_edge] = f;
    ++vertex_visited[face[0]];
    half_edge_mapping[face.back()][face[0]] = next_half_edge;

    const auto starting_half_edge = next_half_edge;
    ++next_half_edge;
    for (UnsignedIndex_t n = 1; n < static_cast<UnsignedIndex_t>(face.size());
         ++n) {
      ending_vertex_mapping[next_half_edge] = face[n];
      previous_half_edge_mapping[next_half_edge] = next_half_edge - 1;
      next_half_edge_mapping[next_half_edge] =
          starting_half_edge +
          (n + 1) % static_cast<UnsignedIndex_t>(face.size());
      face_mapping[next_half_edge] = f;
      half_edge_mapping[face[n - 1]][face[n]] = next_half_edge;
      ++vertex_visited[face[n]];
      ++next_half_edge;
    }
  }

  // Fill in opoosite information for half edges
  next_half_edge = 0;
  for (const auto& face : a_face_brep) {
    opposite_half_edge_mapping[next_half_edge] =
        half_edge_mapping[face[0]][face.back()];
    ++next_half_edge;
    for (UnsignedIndex_t n = 1; n < static_cast<UnsignedIndex_t>(face.size());
         ++n) {
      opposite_half_edge_mapping[next_half_edge] =
          half_edge_mapping[face[n]][face[n - 1]];
      ++next_half_edge;
    }
  }

  const auto max_loc =
      std::max_element(vertex_visited.begin(), vertex_visited.end());
  datum_index = static_cast<UnsignedIndex_t>(
      std::distance(vertex_visited.begin(), max_loc));
  face_triangle_decomposition.reserve(a_face_brep.size() -
                                      static_cast<std::size_t>(*max_loc));

  for (const auto& face : a_face_brep) {
    const auto found_loc = std::find(face.begin(), face.end(), datum_index);
    if (found_loc == face.end()) {
      // Datum not in face, so triangulate and add
      assert(face.size() > 2);
      for (UnsignedIndex_t n = 0; n < face.size() - 2; ++n) {
        face_triangle_decomposition.emplace_back(std::array<UnsignedIndex_t, 3>{
            {face[0], face[n + 1], face[n + 2]}});
      }
    }
  }
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYHEDRONS_POLYHEDRON_CONNECTIVITY_TPP_
