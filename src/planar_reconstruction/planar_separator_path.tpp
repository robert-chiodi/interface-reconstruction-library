// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_PATH_TPP_
#define SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_PATH_TPP_

namespace IRL {

template <template <class NodeType> class NetworkLinkType>
PlanarSeparatorLinkBase<NetworkLinkType>::PlanarSeparatorLinkBase(void)
    : separator_m(nullptr) {}

template <template <class NodeType> class NetworkLinkType>
PlanarSeparatorLinkBase<NetworkLinkType>::PlanarSeparatorLinkBase(
    PlanarSeparator* a_separator)
    : separator_m(a_separator) {
  assert(a_separator != nullptr);
  assert(a_separator->getNumberOfPlanes() < 2);
}

template <template <class NodeType> class NetworkLinkType>
PlanarSeparator&
PlanarSeparatorLinkBase<NetworkLinkType>::getCurrentReconstruction(void) {
  assert(separator_m != nullptr);
  assert(separator_m->getNumberOfPlanes() == 1);
  return separator_m->getCurrentReconstruction();
}

template <template <class NodeType> class NetworkLinkType>
const PlanarSeparator&
PlanarSeparatorLinkBase<NetworkLinkType>::getCurrentReconstruction(void) const {
  assert(separator_m != nullptr);
  assert(separator_m->getNumberOfPlanes() == 1);
  return separator_m->getCurrentReconstruction();
}

template <template <class NodeType> class NetworkLinkType>
inline constexpr NullReconstruction
PlanarSeparatorLinkBase<NetworkLinkType>::getNextReconstruction(void) {
  return NullReconstruction();
}

template <template <class NodeType> class NetworkLinkType>
const PlanarSeparatorLinkBase<NetworkLinkType>*
PlanarSeparatorLinkBase<NetworkLinkType>::getLinkingReconstructionAddress(
    void) const {
  return this;
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_PATH_TPP_
