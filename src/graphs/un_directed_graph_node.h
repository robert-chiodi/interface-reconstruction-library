// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GRAPHS_UN_DIRECTED_GRAPH_NODE_H_
#define SRC_GRAPHS_UN_DIRECTED_GRAPH_NODE_H_

#include <cassert>
#include <vector>

#include "src/parameters/defined_types.h"

namespace IRL {

/// \brief Class to inherit from to provide linking between planes in a
/// planar reconstruction to other reconstructions, where the linking
/// exists PER PLANE.
template <class NodeType>
class UnDirectedGraphNode {
 public:
  using node_type = NodeType;

  /// \brief Default constructor
  UnDirectedGraphNode(void);

  /// \brief Creates an edge for the supplied `a_edge_index`. If this is greater
  /// than the current number of edges, the edge storage will be expanded with
  /// new elements (not at a_edge_index) given nullptr, indicating a
  /// edge to nowhere.
  void setEdgeConnectivity(const UnsignedIndex_t a_edge_index,
                           const NodeType* a_neighbor_ptr);

  /// \brief Return the number of edges. Note, some might be to nowhere, meaning
  /// nullptr.
  UnsignedIndex_t getNumberOfEdges(void) const;

  void setId(const UnsignedIndex_t a_unique_id);
  UnsignedIndex_t getId(void) const;
  bool isIdSet(void) const;

  bool hasNeighbor(const UnsignedIndex_t a_neighbor_index) const;

  /// \brief Return neighboring LocalizedSeparator.
  const NodeType& getNeighbor(const UnsignedIndex_t a_neighbor_index) const;

  const NodeType* getNeighborAddress(
      const UnsignedIndex_t a_neighbor_index) const;

  const NodeType* getNodeMemoryAddress(void) const;

  /// \brief Default destructor.
  ~UnDirectedGraphNode(void) = default;

 private:
  std::vector<const NodeType*> edge_list_m;
  UnsignedIndex_t id_m;
};

}  // namespace IRL

#include "src/graphs/un_directed_graph_node.tpp"

#endif  // SRC_GRAPHS_UN_DIRECTED_GRAPH_NODE_H_
