// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_RESTRICTED_MASKED_LOCALIZED_SEPARATOR_LINK_TPP_
#define SRC_PLANAR_RECONSTRUCTION_RESTRICTED_MASKED_LOCALIZED_SEPARATOR_LINK_TPP_

namespace IRL {

inline RestrictedMaskedLocalizedSeparatorLink::
    RestrictedMaskedLocalizedSeparatorLink(
        const LocalizedSeparatorLink* a_localized_separator_link,
        const PlanarSeparator* a_separator, const WhiteListType* a_white_list)
    : StolenGraph<RestrictedMaskedLocalizedSeparatorLink,
                  LocalizedSeparatorLink>(a_localized_separator_link),
      separator_m(a_separator),
      white_list_m(a_white_list) {}

inline bool RestrictedMaskedLocalizedSeparatorLink::hasNeighbor(
    const UnsignedIndex_t a_neighbor_index) const {
  assert(white_list_m != nullptr);
  if (this->StolenGraph<
          RestrictedMaskedLocalizedSeparatorLink,
          LocalizedSeparatorLink>::hasNeighbor(a_neighbor_index)) {
    return white_list_m->find(this->getNeighborAddress(a_neighbor_index)) !=
           white_list_m->end();
  } else {
    return false;
  }
}

inline RestrictedMaskedLocalizedSeparatorLink
RestrictedMaskedLocalizedSeparatorLink::getNeighbor(
    const UnsignedIndex_t a_neighbor_index) const {
  const LocalizedSeparatorLink* link =
      &this->getStolenNeighbor(a_neighbor_index);
  return RestrictedMaskedLocalizedSeparatorLink(link, separator_m,
                                                white_list_m);
}

inline const PlanarLocalizer&
RestrictedMaskedLocalizedSeparatorLink::getCurrentReconstruction(void) const {
  return this->getStolenGraphNode().getCurrentReconstruction();
}

inline const PlanarSeparator&
RestrictedMaskedLocalizedSeparatorLink::getNextReconstruction(void) const {
  assert(separator_m != nullptr);
  return *separator_m;
}

inline const LocalizedSeparatorLink*
RestrictedMaskedLocalizedSeparatorLink::getLinkingReconstructionAddress(
    void) const {
  return this->getNodeMemoryAddress();
}

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_RESTRICTED_MASKED_LOCALIZED_SEPARATOR_LINK_TPP_
