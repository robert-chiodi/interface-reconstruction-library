// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_LOCALIZED_SEPARATOR_H_
#define SRC_PLANAR_RECONSTRUCTION_LOCALIZED_SEPARATOR_H_

#include "src/planar_reconstruction/joined_reconstructions.h"
#include "src/planar_reconstruction/planar_localizer.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

using LocalizedSeparator =
    JoinedReconstructions<PlanarLocalizer, PlanarSeparator>;
}

#endif  // SRC_PLANAR_RECONSTRUCTION_LOCALIZED_SEPARATOR_H_
