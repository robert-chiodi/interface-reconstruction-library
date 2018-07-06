// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_MASKED_LOCALIZED_SEPARATOR_LINK_TPP_
#define SRC_PLANAR_RECONSTRUCTION_MASKED_LOCALIZED_SEPARATOR_LINK_TPP_

namespace IRL {

inline MaskedLocalizedSeparatorLink::MaskedLocalizedSeparatorLink(
    const LocalizedSeparatorLink* a_localized_separator_link,
    const PlanarSeparator* a_separator)
    : StolenGraph<MaskedLocalizedSeparatorLink, LocalizedSeparatorLink>(
          a_localized_separator_link),
      separator_m(a_separator) {}

inline MaskedLocalizedSeparatorLink MaskedLocalizedSeparatorLink::getNeighbor(
    const UnsignedIndex_t a_neighbor_index) const {
  assert(separator_m != nullptr);
  const LocalizedSeparatorLink* link =
      &this->getStolenNeighbor(a_neighbor_index);
  return MaskedLocalizedSeparatorLink(link, separator_m);
}

inline const PlanarLocalizer&
MaskedLocalizedSeparatorLink::getCurrentReconstruction(void) const {
  return this->getStolenGraphNode().getCurrentReconstruction();
}

inline const PlanarSeparator&
MaskedLocalizedSeparatorLink::getNextReconstruction(void) const {
  assert(separator_m != nullptr);
  return *separator_m;
}

inline const LocalizedSeparatorLink*
MaskedLocalizedSeparatorLink::getLinkingReconstructionAddress(void) const {
  return this->getNodeMemoryAddress();
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_MASKED_LOCALIZED_SEPARATOR_LINK_TPP_
