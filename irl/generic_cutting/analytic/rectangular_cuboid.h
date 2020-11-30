// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_ANALYTIC_RECTANGULAR_CUBOID_H_
#define IRL_GENERIC_CUTTING_ANALYTIC_RECTANGULAR_CUBOID_H_

#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/moments/volume.h"

namespace IRL {

Volume getAnalyticVolume(const RectangularCuboid& a_geometry,
                         const Plane& a_plane);

}  // namespace IRL

#endif // IRL_GENERIC_CUTTING_ANALYTIC_RECTANGULAR_CUBOID_H_
