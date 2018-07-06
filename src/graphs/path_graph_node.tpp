// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GRAPHS_PATH_GRAPH_NODE_TPP_
#define SRC_GRAPHS_PATH_GRAPH_NODE_TPP_

namespace IRL {

template <class NodeType>
PathGraphNode<NodeType>::PathGraphNode(void)
    : next_neighbor_m(nullptr), id_m(static_cast<UnsignedIndex_t>(-1)) {}

template <class NodeType>
void PathGraphNode<NodeType>::setEdgeConnectivity(
    const NodeType* a_neighbor_ptr) {
  next_neighbor_m = a_neighbor_ptr;
}

template <class NodeType>
void PathGraphNode<NodeType>::setId(const UnsignedIndex_t a_id) {
  id_m = a_id;
}

template <class NodeType>
UnsignedIndex_t PathGraphNode<NodeType>::getId(void) const {
  return id_m;
}

template <class NodeType>
bool PathGraphNode<NodeType>::isIdSet(void) const {
  return id_m != static_cast<UnsignedIndex_t>(-1);
}

template <class NodeType>
bool PathGraphNode<NodeType>::hasNeighbor(
    const UnsignedIndex_t a_neighbor_index) const {
  return next_neighbor_m != nullptr;
}

template <class NodeType>
const NodeType& PathGraphNode<NodeType>::getNeighbor(
    const UnsignedIndex_t a_neighbor_index) const {
  assert(next_neighbor_m != nullptr);
  return *(next_neighbor_m);
}

template <class NodeType>
const NodeType* PathGraphNode<NodeType>::getNeighborAddress(
    const UnsignedIndex_t a_neighbor_index) const {
  return next_neighbor_m;
}

template <class NodeType>
const NodeType* PathGraphNode<NodeType>::getNodeMemoryAddress(void) const {
  return static_cast<const NodeType*>(this);
}

}  // namespace IRL

#endif  // SRC_GRAPHS_PATH_GRAPH_NODE_TPP_
