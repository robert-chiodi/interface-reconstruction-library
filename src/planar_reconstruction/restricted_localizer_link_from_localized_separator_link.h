// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_RESTRICTED_LOCALIZER_LINK_FROM_LOCALIZED_SEPARATOR_LINK_H_
#define SRC_PLANAR_RECONSTRUCTION_RESTRICTED_LOCALIZER_LINK_FROM_LOCALIZED_SEPARATOR_LINK_H_

#include <unordered_map>

#include "src/graphs/stolen_graph.h"
#include "src/planar_reconstruction/planar_localizer.h"
#include "src/planar_reconstruction/null_reconstruction.h"
#include "src/planar_reconstruction/localized_separator_link.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

class RestrictedLocalizerLinkFromLocalizedSeparatorLink
    : public StolenGraph<MaskedLocalizedSeparatorLink, LocalizedSeparatorLink> {
  using WhiteListType = std::unordered_map<const LocalizedSeparatorLink*, bool>;

 public:
  /// \brief Default constructor
  RestrictedLocalizerLinkFromLocalizedSeparatorLink(void) = delete;

  /// \brief Construct with pointer to LocalizedSeparatorLink.
  // Also have pointer to a_separator that will always be used
  // to divide the localized regions
  RestrictedLocalizerLinkFromLocalizedSeparatorLink(
      const LocalizedSeparatorLink* a_localized_separator_link,
      const WhiteListType* a_white_list);

  bool hasNeighbor(const UnsignedIndex_t a_neighbor_index) const;

  RestrictedLocalizerLinkFromLocalizedSeparatorLink getNeighbor(
      const UnsignedIndex_t a_neighbor_index) const;

  const PlanarLocalizer& getCurrentReconstruction(void) const;

  static constexpr NullReconstruction getNextReconstruction(void);

  /// \brief Return address of the thing that is linked.
  const LocalizedSeparatorLink* getLinkingReconstructionAddress(void) const;

  /// \brief Default destructor.
  ~RestrictedLocalizerLinkFromLocalizedSeparatorLink(void) = default;

 private:
  const WhiteListType* white_list_m;
};

}  // namespace IRL

#include "src/planar_reconstruction/restricted_localizer_link_from_localized_separator_link.tpp"

#endif  // SRC_PLANAR_RECONSTRUCTION_RESTRICTED_LOCALIZER_LINK_FROM_LOCALIZED_SEPARATOR_LINK_H_
