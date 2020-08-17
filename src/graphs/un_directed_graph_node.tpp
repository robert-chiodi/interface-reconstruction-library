// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GRAPHS_UN_DIRECTED_GRAPH_NODE_TPP_
#define SRC_GRAPHS_UN_DIRECTED_GRAPH_NODE_TPP_

namespace IRL {

template <class NodeType>
UnDirectedGraphNode<NodeType>::UnDirectedGraphNode(void)
    : id_m(static_cast<UnsignedIndex_t>(-1)) {}

template <class NodeType>
void UnDirectedGraphNode<NodeType>::setEdgeConnectivity(
    const UnsignedIndex_t a_edge_index, const NodeType* a_neighbor_ptr) {
  if (a_edge_index >= this->getNumberOfEdges()) {
    edge_list_m.resize(a_edge_index + 1, nullptr);
  }
  edge_list_m[a_edge_index] = a_neighbor_ptr;
}

template <class NodeType>
UnsignedIndex_t UnDirectedGraphNode<NodeType>::getNumberOfEdges(void) const {
  return static_cast<UnsignedIndex_t>(edge_list_m.size());
}

template <class NodeType>
void UnDirectedGraphNode<NodeType>::setId(const UnsignedIndex_t a_unique_id) {
  id_m = a_unique_id;
}

template <class NodeType>
UnsignedIndex_t UnDirectedGraphNode<NodeType>::getId(void) const {
  return id_m;
}

template <class NodeType>
bool UnDirectedGraphNode<NodeType>::isIdSet(void) const {
  return id_m != static_cast<UnsignedIndex_t>(-1);
}

template <class NodeType>
bool UnDirectedGraphNode<NodeType>::hasNeighbor(
    const UnsignedIndex_t a_neighbor_index) const {
  return a_neighbor_index < this->getNumberOfEdges() &&
         edge_list_m[a_neighbor_index] != nullptr;
}

template <class NodeType>
const NodeType& UnDirectedGraphNode<NodeType>::getNeighbor(
    const UnsignedIndex_t a_neighbor_index) const {
  assert(a_neighbor_index < edge_list_m.size());
  assert(edge_list_m[a_neighbor_index] != nullptr);
  return *(edge_list_m[a_neighbor_index]);
}

template <class NodeType>
const NodeType* UnDirectedGraphNode<NodeType>::getNeighborAddress(
    const UnsignedIndex_t a_neighbor_index) const {
  return edge_list_m[a_neighbor_index];
}

template <class NodeType>
const NodeType* UnDirectedGraphNode<NodeType>::getNodeMemoryAddress(
    void) const {
  return static_cast<const NodeType*>(this);
}

}  // namespace IRL

#endif  // SRC_GRAPHS_UN_DIRECTED_GRAPH_NODE_TPP_
