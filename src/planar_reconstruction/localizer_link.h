// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_LOCALIZER_LINK_H_
#define SRC_PLANAR_RECONSTRUCTION_LOCALIZER_LINK_H_

#include "src/planar_reconstruction/localized_separator.h"
#include "src/planar_reconstruction/planar_localizer.h"
#include "src/planar_reconstruction/null_reconstruction.h"
#include "src/planar_reconstruction/localized_separator_link.h"
#include "src/graphs/un_directed_graph_node.h"

namespace IRL {

using LocalizerLink = ReconstructionLink<
    JoinedReconstructions<PlanarLocalizer, NullReconstruction>,
    UnDirectedGraphNode>;

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_LOCALIZER_LINK_H_
