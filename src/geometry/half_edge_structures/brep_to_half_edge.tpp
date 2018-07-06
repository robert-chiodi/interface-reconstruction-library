// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_HALF_EDGE_STRUCTURES_BREP_TO_HALF_EDGE_TPP_
#define SRC_GEOMETRY_HALF_EDGE_STRUCTURES_BREP_TO_HALF_EDGE_TPP_

namespace IRL {

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
template <class PtContainerType, class FaceBREPType>
HalfEdgePolyhedron<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                   kMaxVertices, kMaxFaces>
BREPToHalfEdge<
    PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges, kMaxVertices,
    kMaxFaces>::generateHalfEdgeVersion(const PtContainerType& a_pt_container,
                                        const FaceBREPType& a_face_brep) {
  HalfEdgePolyhedron<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                     kMaxVertices, kMaxFaces>
      half_edge_version;
  setHalfEdgeVersion(a_pt_container, a_face_brep, &half_edge_version);
  return half_edge_version;
}

namespace BREPToHalfEdge_details {

template <class IterType>
IterType checkEdgeExistence(IterType a_begin, IterType a_end,
                            const UnsignedIndex_t a_index) {
  auto count = 0;
  while (a_begin != a_end) {
    if ((*a_begin).first == a_index) {
      return a_begin;
    }
    ++a_begin;
    ++count;
  }
  return a_end;
}
}  // namespace BREPToHalfEdge_details

template <class PtType, class VertexType, class HalfEdgeType, class FaceType,
          UnsignedIndex_t kMaxHalfEdges, UnsignedIndex_t kMaxVertices,
          UnsignedIndex_t kMaxFaces>
template <class PtContainerType, class FaceBREPType,
          class HalfEdgePolyhedronType>
void BREPToHalfEdge<PtType, VertexType, HalfEdgeType, FaceType, kMaxHalfEdges,
                    kMaxVertices, kMaxFaces>::
    setHalfEdgeVersion(const PtContainerType& a_pt_container,
                       const FaceBREPType& a_face_brep,
                       HalfEdgePolyhedronType* a_half_edge_version) {
  // Note: Could also be an unordered_map. Would be better for space, worse for
  // speed.
  std::vector<std::vector<std::pair<UnsignedIndex_t, HalfEdgeType*>>>
      half_edge_mapping(a_pt_container.size());
  for (auto& member : half_edge_mapping) {
    member.reserve(4);
  }

  // E = F + V - 2
  const auto number_of_edges = static_cast<UnsignedIndex_t>(
      a_face_brep.size() + a_pt_container.size() - 2);
  a_half_edge_version->resize(
      2 * number_of_edges, static_cast<UnsignedIndex_t>(a_pt_container.size()),
      static_cast<UnsignedIndex_t>(a_face_brep.size()));

  for (UnsignedIndex_t v = 0;
       v < static_cast<UnsignedIndex_t>(a_pt_container.size()); ++v) {
    a_half_edge_version->getVertex(v).setLocation(a_pt_container[v]);
  }

  UnsignedIndex_t next_half_edge_number = 0;
  for (UnsignedIndex_t f = 0;
       f < static_cast<UnsignedIndex_t>(a_face_brep.size()); ++f) {
    const auto& face = a_face_brep[f];
    auto face_adding = &a_half_edge_version->getFace(f);
    assert(face.size() >= 3);

    // First half edge
    auto current_half_edge =
        &a_half_edge_version->getHalfEdge(next_half_edge_number);

    *current_half_edge = HalfEdgeType(
        &a_half_edge_version->getVertex(face[0]),
        &a_half_edge_version->getHalfEdge(
            next_half_edge_number + static_cast<UnsignedIndex_t>(face.size()) -
            1),
        &a_half_edge_version->getHalfEdge(next_half_edge_number + 1),
        face_adding);

    auto iter = BREPToHalfEdge_details::checkEdgeExistence(
        half_edge_mapping[face[0]].begin(), half_edge_mapping[face[0]].end(),
        face.back());

    if (iter != half_edge_mapping[face[0]].end()) {
      current_half_edge->setOppositeHalfEdge((*iter).second);
      (*iter).second->setOppositeHalfEdge(current_half_edge);
    } else {
      half_edge_mapping[face.back()].emplace_back(
          std::make_pair(face[0], current_half_edge));
    }
    current_half_edge->getVertex()->setHalfEdge(current_half_edge);
    current_half_edge->getFace()->setStartingHalfEdge(current_half_edge);

    // Middle half edges
    decltype(current_half_edge) previous_half_edge;
    decltype(current_half_edge) next_half_edge;
    for (UnsignedIndex_t n = 1;
         n < static_cast<UnsignedIndex_t>(face.size() - 1); ++n) {
      previous_half_edge = current_half_edge;
      next_half_edge =
          &a_half_edge_version->getHalfEdge(next_half_edge_number + n + 1);

      current_half_edge = current_half_edge->getNextHalfEdge();
      *current_half_edge =
          HalfEdgeType(&a_half_edge_version->getVertex(face[n]),
                       previous_half_edge, next_half_edge, face_adding);
      current_half_edge->getVertex()->setHalfEdge(current_half_edge);

      iter = BREPToHalfEdge_details::checkEdgeExistence(
          half_edge_mapping[face[n]].begin(), half_edge_mapping[face[n]].end(),
          face[n - 1]);

      if (iter != half_edge_mapping[face[n]].end()) {
        current_half_edge->setOppositeHalfEdge((*iter).second);
        (*iter).second->setOppositeHalfEdge(current_half_edge);
      } else {
        half_edge_mapping[face[n - 1]].emplace_back(
            std::make_pair(face[n], current_half_edge));
      }
    }

    // Last Half Edges
    previous_half_edge = current_half_edge;
    next_half_edge = &a_half_edge_version->getHalfEdge(next_half_edge_number);
    current_half_edge = current_half_edge->getNextHalfEdge();
    *current_half_edge =
        HalfEdgeType(&a_half_edge_version->getVertex(face.back()),
                     previous_half_edge, next_half_edge, face_adding);
    current_half_edge->getVertex()->setHalfEdge(current_half_edge);

    iter = BREPToHalfEdge_details::checkEdgeExistence(
        half_edge_mapping[face.back()].begin(),
        half_edge_mapping[face.back()].end(), face[face.size() - 2]);

    if (iter != half_edge_mapping[face.back()].end()) {
      current_half_edge->setOppositeHalfEdge((*iter).second);
      (*iter).second->setOppositeHalfEdge(current_half_edge);
    } else {
      half_edge_mapping[face[face.size() - 2]].emplace_back(
          std::make_pair(face.back(), current_half_edge));
    }

    // Increment to next starting half edge number
    next_half_edge_number += static_cast<UnsignedIndex_t>(face.size());
  }
}
}  // namespace IRL

#endif  // SRC_GEOMETRY_HALF_EDGE_STRUCTURES_BREP_TO_HALF_EDGE_TPP_
