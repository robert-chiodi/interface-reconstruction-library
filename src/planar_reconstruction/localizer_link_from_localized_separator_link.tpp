// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry
// operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_LOCALIZER_LINK_FROM_LOCALIZED_SEPARATOR_LINK_TPP_
#define SRC_PLANAR_RECONSTRUCTION_LOCALIZER_LINK_FROM_LOCALIZED_SEPARATOR_LINK_TPP_

namespace IRL {

inline LocalizerLinkFromLocalizedSeparatorLink::
    LocalizerLinkFromLocalizedSeparatorLink(
        const LocalizedSeparatorLink* a_localized_separator_link)
    : StolenGraph<LocalizerLinkFromLocalizedSeparatorLink,
                  LocalizedSeparatorLink>(a_localized_separator_link) {}

inline LocalizerLinkFromLocalizedSeparatorLink
LocalizerLinkFromLocalizedSeparatorLink::getNeighbor(
    const UnsignedIndex_t a_neighbor_index) const {
  const LocalizedSeparatorLink* link =
      &this->getStolenNeighbor(a_neighbor_index);
  return LocalizerLinkFromLocalizedSeparatorLink(link);
}

inline const PlanarLocalizer&
LocalizerLinkFromLocalizedSeparatorLink::getCurrentReconstruction(void) const {
  return this->getStolenGraphNode().getCurrentReconstruction();
}

inline constexpr NullReconstruction
LocalizerLinkFromLocalizedSeparatorLink::getNextReconstruction(void) {
  return NullReconstruction();
}

inline const LocalizedSeparatorLink*
LocalizerLinkFromLocalizedSeparatorLink::getLinkingReconstructionAddress(
    void) const {
  return this->getNodeMemoryAddress();
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_LOCALIZER_LINK_FROM_LOCALIZED_SEPARATOR_LINK_TPP_
