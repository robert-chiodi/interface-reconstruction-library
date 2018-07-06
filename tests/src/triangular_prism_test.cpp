// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/polyhedrons/triangular_prism.h"

#include "gtest/gtest.h"

#include <algorithm>

#include "src/geometry/general/pt.h"

namespace {

using namespace IRL;

TEST(Polyhedra, TriangularPrism) {
  Pt tri_prism_vertices[6] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 1.0, 0.0),
                           Pt(0.0, 0.0, 1.0), Pt(-1.0,0.0,0.0),
						   Pt(-1.0, 1.0, 0.0),Pt(-1.0,0.0,1.0)};
  auto tri_prism = TriangularPrism::fromRawPtPointer(6, tri_prism_vertices);
  Pt centroid = tri_prism.calculateCentroid();

  double correct_volume = 0.5;
  EXPECT_NEAR(tri_prism.calculateVolume(), correct_volume, 1.0e-15);
  EXPECT_NEAR(centroid[0], -0.5, 1.0e-15);
  EXPECT_NEAR(centroid[1], 1.0/3.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], 1.0/3.0, 1.0e-15);

  auto negative_prism = TriangularPrism({Pt(-1.0,0.0,0.0),
	   Pt(-1.0, 1.0, 0.0),Pt(-1.0,0.0,1.0),Pt(0.0, 0.0, 0.0), Pt(0.0, 1.0, 0.0),
       Pt(0.0, 0.0, 1.0)});
  centroid = negative_prism.calculateCentroid();
  correct_volume = -0.5;
  EXPECT_NEAR(negative_prism.calculateVolume(), correct_volume, 1.0e-15);
  EXPECT_NEAR(centroid[0], -0.5, 1.0e-15);
  EXPECT_NEAR(centroid[1], 1.0/3.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], 1.0/3.0, 1.0e-15);
}

}  // namespace
