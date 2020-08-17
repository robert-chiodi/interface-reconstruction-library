// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/interface_reconstruction_methods/reconstruction_cleaning.h"

#include "gtest/gtest.h"
#include "src/generic_cutting/recursive_simplex_cutting/lookup_tables.h"

#include "src/geometry/general/plane.h"
#include "src/helpers/helper.h"
#include "src/helpers/mymath.h"
#include "src/planar_reconstruction/planar_localizer.h"
#include "src/planar_reconstruction/planar_reconstruction.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace {

using namespace IRL;

TEST(ReconstructionCleaning, ReconstructionCleaning) {
  Plane plane_1, plane_2;

  // Ability to clean out bad solutions
  plane_1.normal() = Normal(1.0, 0.0, 0.0);
  plane_1.distance() = (1.0);
  plane_2.normal() = Normal(-1.0, 0.0, 0.0);
  plane_2.distance() = (-1.0);
  PlanarSeparator clean_this =
      PlanarSeparator::fromTwoPlanes(plane_1, plane_2, 1.0);
  clean_this.removePlane(0);
  EXPECT_DOUBLE_EQ(clean_this[0].distance(), -1.0);
  EXPECT_EQ(clean_this.getNumberOfPlanes(), 1);
  clean_this.addPlane(plane_1);
  EXPECT_DOUBLE_EQ(clean_this[1].distance(), 1.0);
  EXPECT_EQ(clean_this.getNumberOfPlanes(), 2);
  clean_this.removePlane(1);
  EXPECT_DOUBLE_EQ(clean_this[0].distance(), -1.0);
  EXPECT_EQ(clean_this.getNumberOfPlanes(), 1);
  clean_this.removePlane(0);
  EXPECT_EQ(clean_this.getNumberOfPlanes(), 0);

  plane_1.distance() = (0.25);
  plane_2.distance() = (-0.25);
  clean_this.addPlane(plane_1);
  clean_this.addPlane(plane_2);
  cleanReconstruction(unit_cell, 0.5, &clean_this);
  EXPECT_EQ(clean_this.getNumberOfPlanes(), 2);

  plane_1.distance() = (10.0);
  clean_this[0] = plane_1;
  cleanReconstruction(unit_cell, 0.25, &clean_this);
  EXPECT_EQ(clean_this.getNumberOfPlanes(), 1);
  EXPECT_DOUBLE_EQ(clean_this[0].distance(), -0.25);
  clean_this.addPlane(plane_1);
  cleanReconstruction(unit_cell, 0.25, &clean_this);
  EXPECT_EQ(clean_this.getNumberOfPlanes(), 1);
  EXPECT_DOUBLE_EQ(clean_this[0].distance(), -0.25);
  clean_this[0] = plane_1;
  cleanReconstruction(unit_cell, 0.5, &clean_this);
  EXPECT_EQ(clean_this.getNumberOfPlanes(), 1);

  plane_1.normal() = Normal(1.0, 0.0, 0.0);
  plane_1.distance() = (0.4);
  plane_2.normal() = Normal(1.0, 0.0, 0.0);
  plane_2.distance() = (0.4);
  PlanarSeparator same_normal =
      PlanarSeparator::fromTwoPlanes(plane_1, plane_2, 1.0);
  cleanReconstruction(unit_cell, 0.75, &same_normal);
  EXPECT_EQ(same_normal.getNumberOfPlanes(), 1);
  EXPECT_DOUBLE_EQ(same_normal[0].distance(), 0.25);

  double cut_volume = getVolumeFraction(unit_cell, same_normal);
  EXPECT_DOUBLE_EQ(cut_volume, 0.75);

  PlanarSeparator flip_setting =
      PlanarSeparator::fromTwoPlanes(Plane(Normal(1.0, 0.0, 0.0), 0.25),
                                     Plane(Normal(0.0, 1.0, 0.0), 0.25), 1.0);
  flip_setting.setFlip(-1.0);
  EXPECT_EQ(flip_setting.isFlipped(), true);

  // Removing out of plane reconstruction for single plane reconstruction
  PlanarSeparator out_of_cell =
      PlanarSeparator::fromOnePlane(Plane(Normal(1.0, 0.0, 0.0), 100.0));
  cleanReconstructionOutOfCell(unit_cell, 0.0, &out_of_cell);
  EXPECT_EQ(out_of_cell.getNumberOfPlanes(), 1);
  out_of_cell.addPlane(Plane(Normal(0.0, 1.0, 0.0), -10000.0));
  cleanReconstruction(unit_cell, 0.0, &out_of_cell);
  EXPECT_EQ(out_of_cell.getNumberOfPlanes(), 1);
  out_of_cell.setNumberOfPlanes(2);
  out_of_cell[0] = (Plane(Normal(0.0, 1.0, 0.0), -10000.0));
  out_of_cell[1] = (Plane(Normal(0.0, 1.0, 0.0), -10000.0));
  cleanReconstructionSameNormal(unit_cell, 0.0, &out_of_cell);
  EXPECT_EQ(out_of_cell.getNumberOfPlanes(), 1);  // Should do nothing for <1>
}
}  // namespace
