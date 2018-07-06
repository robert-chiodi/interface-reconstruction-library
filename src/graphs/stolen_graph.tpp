// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GRAPHS_STOLEN_GRAPH_TPP_
#define SRC_GRAPHS_STOLEN_GRAPH_TPP_

namespace IRL {

template <class NodeType, class StolenGraphType>
StolenGraph<NodeType, StolenGraphType>::StolenGraph(
    const StolenGraphType* a_stolen_graph)
    : stolen_graph_m(a_stolen_graph) {}

template <class NodeType, class StolenGraphType>
UnsignedIndex_t StolenGraph<NodeType, StolenGraphType>::getNumberOfEdges(
    void) const {
  return stolen_graph_m->getNumberOfEdges();
}

template <class NodeType, class StolenGraphType>
UnsignedIndex_t StolenGraph<NodeType, StolenGraphType>::getId(void) const {
  return stolen_graph_m->getId();
}

template <class NodeType, class StolenGraphType>
bool StolenGraph<NodeType, StolenGraphType>::isIdSet(void) const {
  return stolen_graph_m->isIdSet();
}

template <class NodeType, class StolenGraphType>
bool StolenGraph<NodeType, StolenGraphType>::hasNeighbor(
    const UnsignedIndex_t a_neighbor_index) const {
  return stolen_graph_m->hasNeighbor(a_neighbor_index);
}

template <class NodeType, class StolenGraphType>
const typename StolenGraph<NodeType, StolenGraphType>::stolen_node_type&
StolenGraph<NodeType, StolenGraphType>::getStolenNeighbor(
    const UnsignedIndex_t a_neighbor_index) const {
  return stolen_graph_m->getNeighbor(a_neighbor_index);
}

template <class NodeType, class StolenGraphType>
const typename StolenGraph<NodeType, StolenGraphType>::stolen_node_type*
StolenGraph<NodeType, StolenGraphType>::getNeighborAddress(
    const UnsignedIndex_t a_neighbor_index) const {
  return stolen_graph_m->getNeighborAddress(a_neighbor_index);
}

template <class NodeType, class StolenGraphType>
const typename StolenGraph<NodeType, StolenGraphType>::stolen_node_type*
StolenGraph<NodeType, StolenGraphType>::getNodeMemoryAddress(void) const {
  return stolen_graph_m->getNodeMemoryAddress();
}

template <class NodeType, class StolenGraphType>
const typename StolenGraph<NodeType, StolenGraphType>::stolen_node_type&
StolenGraph<NodeType, StolenGraphType>::getStolenGraphNode(void) const {
  assert(stolen_graph_m != nullptr);
  return *stolen_graph_m;
}

}  // namespace IRL

#endif  // SRC_GRAPHS_STOLEN_GRAPH_TPP_
