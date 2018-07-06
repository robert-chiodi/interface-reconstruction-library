// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_PATH_H_
#define SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_PATH_H_

#include "src/graphs/path_graph_node.h"
#include "src/planar_reconstruction/localized_separator.h"
#include "src/planar_reconstruction/localized_separator_link.h"
#include "src/planar_reconstruction/null_reconstruction.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

/// \brief This class handles the ability to chain together PlanarSeparator
/// objects. Currently, due to restrictions in the implementation of the cutting
/// routines, these PlanarSeparator objects are limited to a single plane. In
/// actuality, the limitation is on requring that the planes represent a convex
/// volume (meaning flip_cut_m = 1.0). This is not a fundamental limitation, but
/// simply requires some additional work in the cutting routines to enable the
/// use of PlanarSeparators representing non-convex volumes. The necessity to
/// have a single plane is enforced through assertions. Make sure to not run
/// without assertions off and multiple-plane PlanarSeparators!

// NOTE: This could be the definition of PlanarSeparatorPath. For now, will use
// the custom class below it because I can assert single-plane PlanarSeparators.
//using PlanarSeparatorPath = ReconstructionLink<
//    JoinedReconstructions<PlanarSeparator, NullReconstruction>, PathGraphNode>;

 template <template <typename NodeType> class NetworkLinkType>
 class PlanarSeparatorLinkBase
    : public NetworkLinkType<PlanarSeparatorLinkBase<NetworkLinkType>> {
 public:
  /// \brief Default constructor
  PlanarSeparatorLinkBase(void);

  /// \brief Construct with pointer to localized separator. Need to set unique
  /// ID later.
  explicit PlanarSeparatorLinkBase(PlanarSeparator* a_separator);

  PlanarSeparator& getCurrentReconstruction(void);
  const PlanarSeparator& getCurrentReconstruction(void) const;

  static constexpr NullReconstruction getNextReconstruction(void);

  const PlanarSeparatorLinkBase* getLinkingReconstructionAddress(void) const;

  /// \brief Default destructor.
  ~PlanarSeparatorLinkBase(void) = default;

 private:
  PlanarSeparator* separator_m;
};

 using PlanarSeparatorPath = PlanarSeparatorLinkBase<PathGraphNode>;

}  // namespace IRL

#include "src/planar_reconstruction/planar_separator_path.tpp"

#endif  // SRC_PLANAR_RECONSTRUCTION_PLANAR_SEPARATOR_PATH_H_
