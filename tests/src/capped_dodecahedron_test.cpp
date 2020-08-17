// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/polyhedrons/capped_dodecahedron.h"

#include "gtest/gtest.h"

#include <algorithm>

#include "src/geometry/general/pt.h"

namespace {

using namespace IRL;

TEST(Polyhedra, CappedDodecahedron) {
  // Test self intersecting hexahedron

  // In X
  CappedDodecahedron capped_dodeca(
      {Pt(0.5, -0.5, -0.5), Pt(0.5, 0.5, -0.5), Pt(0.5, 0.5, 0.5),
       Pt(0.5, -0.5, 0.5), Pt(-0.5, -0.5, -0.5), Pt(-0.5, 0.5, -0.5),
       Pt(-0.5, 0.5, 0.5), Pt(-0.5, -0.5, 0.5), Pt(-0.5, 0.0, 0.0)});
  Pt centroid = capped_dodeca.calculateCentroid();
  EXPECT_NEAR(capped_dodeca.calculateVolume(), 1.0, 1.0e-14);
  EXPECT_NEAR(centroid[0], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[1], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], 0.0, 1.0e-15);

  capped_dodeca.adjustCapToMatchVolume(2.0);
  centroid = capped_dodeca.calculateCentroid();
  EXPECT_NEAR(capped_dodeca.calculateVolume(), 2.0, 1.0e-14);
  EXPECT_NEAR(capped_dodeca[8][0], -3.5, 1.0e-14);
  EXPECT_NEAR(capped_dodeca[8][1], 0.0, 1.0e-14);
  EXPECT_NEAR(capped_dodeca[8][2], 0.0, 1.0e-14);
  EXPECT_NEAR(centroid[0], -5.0 / 8.0, 1.0e-15);
  EXPECT_NEAR(centroid[1], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], 0.0, 1.0e-15);

  // In -X
  capped_dodeca = CappedDodecahedron(
      {Pt(-0.5, -0.5, -0.5), Pt(-0.5, 0.5, -0.5), Pt(-0.5, 0.5, 0.5),
       Pt(-0.5, -0.5, 0.5), Pt(0.5, -0.5, -0.5), Pt(0.5, 0.5, -0.5),
       Pt(0.5, 0.5, 0.5), Pt(0.5, -0.5, 0.5), Pt(0.5, 0.0, 0.0)});
  centroid = capped_dodeca.calculateCentroid();
  EXPECT_NEAR(capped_dodeca.calculateVolume(), -1.0, 1.0e-14);
  EXPECT_NEAR(centroid[0], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[1], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], 0.0, 1.0e-15);

  capped_dodeca.adjustCapToMatchVolume(-2.0);
  centroid = capped_dodeca.calculateCentroid();
  EXPECT_NEAR(capped_dodeca.calculateVolume(), -2.0, 1.0e-14);
  EXPECT_NEAR(capped_dodeca[8][0], 3.5, 1.0e-14);
  EXPECT_NEAR(capped_dodeca[8][1], 0.0, 1.0e-14);
  EXPECT_NEAR(capped_dodeca[8][2], 0.0, 1.0e-14);
  EXPECT_NEAR(centroid[0], 5.0 / 8.0, 1.0e-15);
  EXPECT_NEAR(centroid[1], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], 0.0, 1.0e-15);

  // In Y
  capped_dodeca = CappedDodecahedron(
      {Pt(0.5, 0.5, -0.5), Pt(-0.5, 0.5, -0.5), Pt(-0.5, 0.5, 0.5),
       Pt(0.5, 0.5, 0.5), Pt(0.5, -0.5, -0.5), Pt(-0.5, -0.5, -0.5),
       Pt(-0.5, -0.5, 0.5), Pt(0.5, -0.5, 0.5), Pt(0.0, -0.5, 0.0)});
  centroid = capped_dodeca.calculateCentroid();
  EXPECT_NEAR(capped_dodeca.calculateVolume(), 1.0, 1.0e-14);
  EXPECT_NEAR(centroid[0], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[1], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], 0.0, 1.0e-15);

  capped_dodeca.adjustCapToMatchVolume(2.0);
  centroid = capped_dodeca.calculateCentroid();
  EXPECT_NEAR(capped_dodeca.calculateVolume(), 2.0, 1.0e-14);
  EXPECT_NEAR(capped_dodeca[8][0], 0.0, 1.0e-14);
  EXPECT_NEAR(capped_dodeca[8][1], -3.5, 1.0e-14);
  EXPECT_NEAR(capped_dodeca[8][2], 0.0, 1.0e-14);
  EXPECT_NEAR(centroid[0], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[1], -5.0 / 8.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], 0.0, 1.0e-15);

  // In -Z
  capped_dodeca = CappedDodecahedron(
      {Pt(0.5, 0.5, -0.5), Pt(-0.5, 0.5, -0.5), Pt(-0.5, -0.5, -0.5),
       Pt(0.5, -0.5, -0.5), Pt(0.5, 0.5, 0.5), Pt(-0.5, 0.5, 0.5),
       Pt(-0.5, -0.5, 0.5), Pt(0.5, -0.5, 0.5), Pt(0.0, 0.0, 0.5)});
  centroid = capped_dodeca.calculateCentroid();
  EXPECT_NEAR(capped_dodeca.calculateVolume(), -1.0, 1.0e-14);
  EXPECT_NEAR(centroid[0], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[1], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], 0.0, 1.0e-15);

  capped_dodeca.adjustCapToMatchVolume(2.0);
  centroid = capped_dodeca.calculateCentroid();
  EXPECT_NEAR(capped_dodeca.calculateVolume(), 2.0, 1.0e-14);
  EXPECT_NEAR(capped_dodeca[8][0], 0.0, 1.0e-14);
  EXPECT_NEAR(capped_dodeca[8][1], 0.0, 1.0e-14);
  EXPECT_NEAR(capped_dodeca[8][2], -8.5, 1.0e-14);
  EXPECT_NEAR(centroid[0], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[1], 0.0, 1.0e-15);
  EXPECT_NEAR(centroid[2], -2.625, 1.0e-15);
}

}  // namespace
