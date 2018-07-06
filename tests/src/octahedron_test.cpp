// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/polyhedrons/octahedron.h"

#include "gtest/gtest.h"

#include <algorithm>

#include "src/geometry/general/pt.h"

namespace {

using namespace IRL;

TEST(Polyhedra, Octahedron) {
  // Degenerate case when returns to
  Pt octahedron_vertices[6] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 1.0, 0.0),
                           Pt(0.0, 0.0, 1.0), Pt(-1.0,0.0,0.0),
						   Pt(-1.0, 1.0, 0.0),Pt(-1.0,0.0,1.0)};
  auto octahedron = Octahedron::fromRawPtPointer(6, octahedron_vertices);
  Pt centroid = octahedron.calculateCentroid();

  double correct_volume = 0.5;
  EXPECT_NEAR(octahedron.calculateVolume(), correct_volume, 1.0e-15);
  EXPECT_NEAR(centroid[0], -0.5, 1.0e-15);
  EXPECT_NEAR(centroid[1], 1.0/3.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], 1.0/3.0, 1.0e-15);

  auto negative_oct = Octahedron({Pt(-1.0,0.0,0.0),
	   Pt(-1.0, 1.0, 0.0),Pt(-1.0,0.0,1.0),Pt(0.0, 0.0, 0.0), Pt(0.0, 1.0, 0.0),
       Pt(0.0, 0.0, 1.0)});
  centroid = negative_oct.calculateCentroid();
  correct_volume = -0.5;
  EXPECT_NEAR(negative_oct.calculateVolume(), correct_volume, 1.0e-15);
  EXPECT_NEAR(centroid[0], -0.5, 1.0e-15);
  EXPECT_NEAR(centroid[1], 1.0/3.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], 1.0/3.0, 1.0e-15);
}

}  // namespace
