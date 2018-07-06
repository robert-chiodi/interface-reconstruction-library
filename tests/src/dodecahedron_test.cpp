// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/polyhedrons/dodecahedron.h"

#include "gtest/gtest.h"

#include <algorithm>

#include "src/geometry/general/pt.h"

namespace {

using namespace IRL;

TEST(Polyhedra, Dodecahedron) {
  // Test self intersecting hexahedron
  Pt hex_vertex_list[8] = {Pt(0.0, 0.0, 0.0),  Pt(0.0, 1.0, 0.0),
                           Pt(0.0, 1.0, 1.0),  Pt(0.0, 0.0, 1.0),
                           Pt(-1.0, 0.0, 0.0), Pt(0.5, 1.0, 0.0),
                           Pt(0.5, 1.0, 1.0),  Pt(-1.0, 0.0, 1.0)};
  Dodecahedron dodeca = Dodecahedron::fromRawPtPointer(8, hex_vertex_list);
  Pt centroid = dodeca.calculateCentroid();

  double correct_volume = 0.25;
  EXPECT_NEAR(dodeca.calculateVolume(), correct_volume, 1.0e-15);
  EXPECT_NEAR(centroid[0], -0.5, 1.0e-15);
  EXPECT_NEAR(centroid[1], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], 0.5, 1.0e-15);
}

}  // namespace
