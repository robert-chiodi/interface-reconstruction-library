// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_RECONSTRUCTION_LINK_TPP_
#define SRC_PLANAR_RECONSTRUCTION_RECONSTRUCTION_LINK_TPP_

namespace IRL {

template <class ReconstructionType, template <class NodeType> class GraphType>
const ReconstructionLink<ReconstructionType, GraphType>*
ReconstructionLink<ReconstructionType,
                   GraphType>::getLinkingReconstructionAddress(void) const {
  return this;
}

template <class ReconstructionType, template <class NodeType> class GraphType>
std::ostream& operator<<(
    std::ostream& out,
    const ReconstructionLink<ReconstructionType, GraphType>& a_reconstruction) {
  out << "First reconstruction is " << std::endl;
  out << a_reconstruction.getCurrentReconstruction() << std::endl;
  out << "Second reconstruction is " << std::endl;
  out << a_reconstruction.getNextReconstruction() << std::endl;

  if(a_reconstruction.isIdSet()){
    out << "This link's ID is " << a_reconstruction.getId() << std::endl;
  } else{
	out << "This link has NO ID" << std::endl;
  }


  out << "Connectivity to IDs of " << std::endl;
  for (UnsignedIndex_t plane = 0;
       plane < a_reconstruction.getCurrentReconstruction().getNumberOfPlanes();
       ++plane) {
    if (a_reconstruction.hasNeighbor(plane)) {
      out << "Plane " << plane << " has a neighbor with global ID "
          << a_reconstruction.getNeighbor(plane).getId()
          << " at memory address " << a_reconstruction.getNeighborAddress(plane)
          << std::endl;
    } else {
      out << "Plane " << plane << " has NO neighbor " << std::endl;
    }
  }
  return out;
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_RECONSTRUCTION_LINK_TPP_
