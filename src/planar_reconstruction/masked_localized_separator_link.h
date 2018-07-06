// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_MASKED_LOCALIZED_SEPARATOR_LINK_H_
#define SRC_PLANAR_RECONSTRUCTION_MASKED_LOCALIZED_SEPARATOR_LINK_H_

#include "src/graphs/stolen_graph.h"
#include "src/planar_reconstruction/planar_localizer.h"
#include "src/planar_reconstruction/localized_separator_link.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

class MaskedLocalizedSeparatorLink
    : public StolenGraph<MaskedLocalizedSeparatorLink, LocalizedSeparatorLink> {
 public:
  MaskedLocalizedSeparatorLink(void) = delete;

  /// \brief Construct with pointer to LocalizedSeparatorLink.
  // Also have pointer to a_separator that will always be used
  // to divide the localized regions
  MaskedLocalizedSeparatorLink(
      const LocalizedSeparatorLink* a_localized_separator_link,
      const PlanarSeparator* a_separator);

  MaskedLocalizedSeparatorLink getNeighbor(
      const UnsignedIndex_t a_neighbor_index) const;

  const PlanarLocalizer& getCurrentReconstruction(void) const;

  const PlanarSeparator& getNextReconstruction(void) const;

  /// \brief Return address of the thing that is linked.
  const LocalizedSeparatorLink* getLinkingReconstructionAddress(void) const;

  /// \brief Default destructor.
  ~MaskedLocalizedSeparatorLink(void) = default;

 private:
  const PlanarSeparator* separator_m;
};

}  // namespace IRL

#include "src/planar_reconstruction/masked_localized_separator_link.tpp"

#endif  // SRC_PLANAR_RECONSTRUCTION_MASKED_LOCALIZED_SEPARATOR_LINK_H_
