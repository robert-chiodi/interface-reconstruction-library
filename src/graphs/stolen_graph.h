// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GRAPHS_STOLEN_GRAPH_H_
#define SRC_GRAPHS_STOLEN_GRAPH_H_

#include <cassert>
#include <vector>

#include "src/parameters/defined_types.h"

namespace IRL {

/// \brief Class to inherit from that will utilize an already set up graph
/// and mask some of its behaviors by shadowing in the NodeType class. In order
/// to have a consistent interface with other graph classes, the NodeType class
/// that inherits StolenGraph (which works due to CRTP) must implement a
/// NodeType getNeighbor(const UnsignedIndex_t a_neighbor_index) const method
/// See planar_reconstruction/localizer_link_from_localized_separator_link.h
/// as an example.
template <class NodeType, class StolenGraphType>
class StolenGraph {
  using stolen_node_type = typename StolenGraphType::node_type;

 public:
  /// \brief Default constructor
  explicit StolenGraph(const StolenGraphType* a_stolen_graph);

  /// \brief Return the number of edges. Note, some might be to nowhere, meaning
  /// nullptr.
  UnsignedIndex_t getNumberOfEdges(void) const;

  UnsignedIndex_t getId(void) const;
  bool isIdSet(void) const;

  bool hasNeighbor(const UnsignedIndex_t a_neighbor_index) const;

  /// \brief Return neighboring LocalizedSeparator.
  const stolen_node_type& getStolenNeighbor(
      const UnsignedIndex_t a_neighbor_index) const;

  const stolen_node_type* getNeighborAddress(
      const UnsignedIndex_t a_neighbor_index) const;

  const stolen_node_type* getNodeMemoryAddress(void) const;

  const stolen_node_type& getStolenGraphNode(void) const;

  /// \brief Default destructor.
  ~StolenGraph(void) = default;

 private:
  const StolenGraphType* stolen_graph_m;
};

}  // namespace IRL

#include "src/graphs/stolen_graph.tpp"

#endif  // SRC_GRAPHS_STOLEN_GRAPH_H_
