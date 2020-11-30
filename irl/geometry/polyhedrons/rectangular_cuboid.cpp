// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/geometry/polyhedrons/rectangular_cuboid.h"

#include "irl/geometry/general/pt.h"
#include "irl/geometry/polygons/tri.h"
#include "irl/planar_reconstruction/planar_localizer.h"

namespace IRL {

const RectangularCuboid unit_cell =
    RectangularCuboid::fromBoundingPts(Pt(-0.5, -0.5, -0.5), Pt(0.5, 0.5, 0.5));

}  // namespace IRL
