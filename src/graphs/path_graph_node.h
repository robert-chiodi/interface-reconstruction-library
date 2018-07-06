// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GRAPHS_PATH_GRAPH_NODE_H_
#define SRC_GRAPHS_PATH_GRAPH_NODE_H_

#include <cassert>

#include "src/parameters/defined_types.h"

namespace IRL {

/// \brief Class to inherit from to provide linking between reconstructions,
/// where the linking exists PER RECONSTRUCTION. This also follows a
/// directed path-graph (linear graph), with a starting node and constant
/// direction to the end node. The interface is kept the same
/// as for the UnDirectedGraphNode in order to facilitate
/// reuse in the cutting routines. All neighbor indices
/// will simply point to the next node in the path-graph.
template <class NodeType>
class PathGraphNode {
 public:
  using node_type = NodeType;

  /// \brief Default constructor
  PathGraphNode(void);

  /// \brief Supplied the pointer for the neighbor at the index.
  void setEdgeConnectivity(const NodeType* a_neighbor_ptr);

  /// \brief Set unique Id for this node in the graph.
  void setId(const UnsignedIndex_t a_id);
  UnsignedIndex_t getId(void) const;
  bool isIdSet(void) const;

  // Keep receiving a_neighbor_index in order to keep consistent interface with
  // UnDirectedGraphNode. It is optional though.
  bool hasNeighbor(const UnsignedIndex_t a_neighbor_index = 0) const;

  // Keep receiving a_neighbor_index in order to keep consistent interface with
  // UnDirectedGraphNode. It is optional though.
  const NodeType& getNeighbor(const UnsignedIndex_t a_neighbor_index = 0) const;

  // Keep receiving a_neighbor_index in order to keep consistent interface with
  // UnDirectedGraphNode. It is optional though.
  const NodeType* getNeighborAddress(
      const UnsignedIndex_t a_neighbor_index = 0) const;

  const NodeType* getNodeMemoryAddress(void) const;

  /// \brief Default destructor.
  ~PathGraphNode(void) = default;

 private:
  const NodeType* next_neighbor_m;
  UnsignedIndex_t id_m;
};

}  // namespace IRL

#include "src/graphs/path_graph_node.tpp"

#endif  // SRC_GRAPHS_PATH_GRAPH_NODE_H_
