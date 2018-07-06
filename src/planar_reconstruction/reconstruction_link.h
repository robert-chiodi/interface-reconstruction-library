// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_RECONSTRUCTION_LINK_H_
#define SRC_PLANAR_RECONSTRUCTION_RECONSTRUCTION_LINK_H_

#include <ostream>
#include <vector>

#include "src/graphs/un_directed_graph_node.h"

namespace IRL {

/// \brief This is a class template for linking together reconstructions that
/// are already pointers to something else. This then does not require a
/// different constructor because we are already taking pointers of the
/// underlying reconstructions. For classes that are currently storing their
/// own data, it would probably be best to generate a wrapping pointer class,
/// which can be done with JoinedReconstruction<ClassYouWant,
/// NullReconstruction>.
template <class ReconstructionType, template <class NodeType> class GraphType>
class ReconstructionLink
    : public ReconstructionType,
      public GraphType<ReconstructionLink<ReconstructionType, GraphType>> {
 public:
  using ReconstructionType::ReconstructionType;

  const ReconstructionLink* getLinkingReconstructionAddress(void) const;

 private:
};

template <class ReconstructionType, template <class NodeType> class GraphType>
std::ostream& operator<<(
    std::ostream& out,
    const ReconstructionLink<ReconstructionType, GraphType>& a_reconstruction);

}  // namespace IRL

#include "src/planar_reconstruction/reconstruction_link.tpp"

#endif  // SRC_PLANAR_RECONSTRUCTION_RECONSTRUCTION_LINK_H_
