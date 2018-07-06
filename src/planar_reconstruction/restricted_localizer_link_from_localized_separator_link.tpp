// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_RESTRICTED_LOCALIZER_LINK_FROM_LOCALIZED_SEPARATOR_LINK_TPP_
#define SRC_PLANAR_RECONSTRUCTION_RESTRICTED_LOCALIZER_LINK_FROM_LOCALIZED_SEPARATOR_LINK_TPP_

namespace IRL {

inline RestrictedLocalizerLinkFromLocalizedSeparatorLink::
    RestrictedLocalizerLinkFromLocalizedSeparatorLink(
        const LocalizedSeparatorLink* a_localized_separator_link,
        const WhiteListType* a_white_list)
    : StolenGraph<MaskedLocalizedSeparatorLink, LocalizedSeparatorLink>(
          a_localized_separator_link),
      white_list_m(a_white_list) {}

inline bool RestrictedLocalizerLinkFromLocalizedSeparatorLink::hasNeighbor(
    const UnsignedIndex_t a_neighbor_index) const {
  assert(white_list_m != nullptr);
  if (StolenGraph<MaskedLocalizedSeparatorLink,
                  LocalizedSeparatorLink>::hasNeighbor(a_neighbor_index)) {
    return white_list_m->find(this->getNeighborAddress(a_neighbor_index)) !=
           white_list_m->end();
  } else {
    return false;
  }
}

inline RestrictedLocalizerLinkFromLocalizedSeparatorLink
RestrictedLocalizerLinkFromLocalizedSeparatorLink::getNeighbor(
    const UnsignedIndex_t a_neighbor_index) const {
  const LocalizedSeparatorLink* link =
      &this->getStolenNeighbor(a_neighbor_index);
  return RestrictedLocalizerLinkFromLocalizedSeparatorLink(link, white_list_m);
}

inline const PlanarLocalizer&
RestrictedLocalizerLinkFromLocalizedSeparatorLink::getCurrentReconstruction(
    void) const {
  return this->getStolenGraphNode().getCurrentReconstruction();
}

inline constexpr NullReconstruction
RestrictedLocalizerLinkFromLocalizedSeparatorLink::getNextReconstruction(void) {
  return NullReconstruction();
}

inline const LocalizedSeparatorLink*
RestrictedLocalizerLinkFromLocalizedSeparatorLink::
    getLinkingReconstructionAddress(void) const {
  return this->getNodeMemoryAddress();
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_RESTRICTED_LOCALIZER_LINK_FROM_LOCALIZED_SEPARATOR_LINK_TPP_
