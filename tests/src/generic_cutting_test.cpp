// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "irl/generic_cutting/generic_cutting.h"

#include <vector>

#include "gtest/gtest.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/polygons/divided_polygon.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/geometry/polygons/tri.h"
#include "irl/geometry/polyhedrons/capped_octahedron_variations/capped_octahedron_LLL.h"
#include "irl/helpers/geometric_cutting_helpers.h"
#include "irl/interface_reconstruction_methods/reconstruction_cleaning.h"
#include "irl/moments/accumulate_wrapper.h"
#include "irl/planar_reconstruction/localized_separator.h"
#include "irl/planar_reconstruction/localized_separator_group.h"
#include "irl/planar_reconstruction/localized_separator_group_link.h"
#include "irl/planar_reconstruction/localized_separator_link.h"
#include "irl/planar_reconstruction/localizer_link.h"
#include "irl/planar_reconstruction/planar_separator_path.h"
#include "irl/planar_reconstruction/planar_separator_path_group.h"

#include <iostream>

namespace {

using namespace IRL;

TEST(GenericCutting, cutTetByPlanarReconstruction) {
  Tet tet({Pt(0.0, 1.0, 0.0), Pt(0.0, 0.0, 1.0), Pt(1.0, 0.0, 0.0),
           Pt(0.0, 0.0, -1.0)});
  EXPECT_DOUBLE_EQ(tet.calculateVolume(), 0.5 * 2.0 / 3.0);

  Plane plane_up_y(Normal(0.0, 1.0, 0.0), 0.25);

  Pt midpoint_0 = Pt::fromEdgeIntersection(
      tet[0], plane_up_y.signedDistanceToPoint(tet[0]), tet[1],
      plane_up_y.signedDistanceToPoint(tet[1]));
  Pt midpoint_1 = Pt::fromEdgeIntersection(
      tet[0], plane_up_y.signedDistanceToPoint(tet[0]), tet[2],
      plane_up_y.signedDistanceToPoint(tet[2]));
  Pt midpoint_2 = Pt::fromEdgeIntersection(
      tet[0], plane_up_y.signedDistanceToPoint(tet[0]), tet[3],
      plane_up_y.signedDistanceToPoint(tet[3]));
  Tet upper_tet({tet[0], midpoint_0, midpoint_1, midpoint_2});

  // Case with plane_up_y
  PlanarSeparator recon_one_plane = PlanarSeparator::fromOnePlane(plane_up_y);
  // recon_one_plane.flipCutting();
  SeparatedMoments<VolumeMoments> cut_from_tet =
      getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          tet, recon_one_plane);
  EXPECT_DOUBLE_EQ(tet.calculateVolume() - upper_tet.calculateVolume(),
                   cut_from_tet[0].volume());
  EXPECT_DOUBLE_EQ(upper_tet.calculateVolume(), cut_from_tet[1].volume());

  // Test by using to cut cube by sheet
  // Cut cube by one plane reconstruction
  PlanarSeparator reconstruction_for_cube =
      PlanarSeparator::fromTwoPlanes(Plane(Normal(0.0, 1.0, 0.0), 0.25),
                                     Plane(Normal(0.0, -1.0, 0.0), 0.0), 1.0);

  SeparatedMoments<VolumeMoments> phase_moments;
  for (UnsignedIndex_t t = 0;
       t < unit_cell.getNumberOfSimplicesInDecomposition(); ++t) {
    phase_moments += getVolumeMoments<SeparatedMoments<VolumeMoments>>(
        unit_cell.getSimplexFromDecomposition(t), reconstruction_for_cube);
  }
  double liquid_volume = 0.25;
  phase_moments.normalizeByVolume();
  Pt liquid_centroid(0.0, 0.125, 0.0);
  EXPECT_DOUBLE_EQ(phase_moments[0].volume(), liquid_volume);
  EXPECT_NEAR(phase_moments[0].centroid()[0], liquid_centroid[0], DBL_EPSILON);
  EXPECT_DOUBLE_EQ(phase_moments[0].centroid()[1], liquid_centroid[1]);
  EXPECT_NEAR(phase_moments[0].centroid()[2], liquid_centroid[2], DBL_EPSILON);
}

TEST(GenericCutting, cutRectangularCuboidByPlanes) {
  RectangularCuboid cube = unit_cell;
  Plane plane_up_y(Normal(0.0, 1.0, 0.0), 0.25);
  Plane plane_down_y(Normal(0.0, -1.0, 0.0), 0.25);
  PlanarSeparator recon_one_plane = PlanarSeparator::fromOnePlane(plane_up_y);

  // Cut cube by one plane reconstruction
  SeparatedMoments<VolumeMoments> from_cube_cut =
      getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          cube, recon_one_plane);
  double liquid_volume = 0.75;
  Pt liquid_centroid(0.0, -0.125, 0.0);
  double gas_volume = 0.25;
  Pt gas_centroid(0.0, 0.375, 0.0);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].volume(), liquid_volume);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].centroid()[0], liquid_centroid[0]);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].centroid()[1], liquid_centroid[1]);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].centroid()[2], liquid_centroid[2]);
  EXPECT_DOUBLE_EQ(from_cube_cut[1].volume(), gas_volume);
  EXPECT_DOUBLE_EQ(from_cube_cut[1].centroid()[0], gas_centroid[0]);
  EXPECT_DOUBLE_EQ(from_cube_cut[1].centroid()[1], gas_centroid[1]);
  EXPECT_DOUBLE_EQ(from_cube_cut[1].centroid()[2], gas_centroid[2]);

  // Cut cube by two plane reconstruction
  liquid_volume = 0.5;
  liquid_centroid = Pt(0.0, 0.0, 0.0);
  gas_volume = 0.5;
  gas_centroid = Pt(0.0, 0.0, 0.0);
  PlanarSeparator recon_two_plane =
      PlanarSeparator::fromTwoPlanes(plane_up_y, plane_down_y, 1.0);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, recon_two_plane);
  EXPECT_NEAR(from_cube_cut[0].volume(), liquid_volume, 1.0e-14);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[0].centroid()[2], liquid_centroid[2], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].volume(), gas_volume, 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], 1.0e-14);

  // Cut on shifted cube
  double shift_x = -4.8;
  double shift_y = 0.125;
  double shift_z = 3.28;
  cube.shift(shift_x, shift_y, shift_z);
  liquid_volume = 0.5;
  liquid_centroid = Pt(shift_x, 0.0, shift_z);
  gas_volume = 0.5;
  gas_centroid = Pt(shift_x, 2 * shift_y, shift_z);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, recon_two_plane);
  EXPECT_NEAR(from_cube_cut[0].volume(), liquid_volume, 1.0e-14);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[0].centroid()[2], liquid_centroid[2], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].volume(), gas_volume, 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], 1.0e-14);

  cube.shift(-shift_x, -shift_y, -shift_z);
  Plane plane_up_x(Normal(1.0, 0.0, 0.0), 0.25);
  PlanarSeparator recon_one_plane_right =
      PlanarSeparator::fromOnePlane(plane_up_x);
  liquid_volume = 0.75;
  liquid_centroid = Pt(-0.125, 0.0, 0.0);
  gas_volume = 0.25;
  gas_centroid = Pt(0.375, 0.0, 0.0);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, recon_one_plane_right);
  EXPECT_NEAR(from_cube_cut[0].volume(), liquid_volume, 1.0e-14);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[0].centroid()[2], liquid_centroid[2], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].volume(), gas_volume, 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], 1.0e-14);

  Plane plane_45(Normal::normalized(1.0, 1.0, 0.0), 0.0);
  PlanarSeparator recon_45 = PlanarSeparator::fromOnePlane(plane_45);
  liquid_volume = 0.5;
  liquid_centroid = Pt(-0.5 / 3.0, -0.5 / 3.0, 0.0);
  gas_volume = 0.5;
  gas_centroid = Pt(0.5 / 3.0, 0.5 / 3.0, 0.0);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, recon_45);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].volume(), liquid_volume);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[2], liquid_centroid[2], DBL_EPSILON);
  EXPECT_DOUBLE_EQ(from_cube_cut[1].volume(), gas_volume);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], DBL_EPSILON);
}

TEST(GenericCutting, getRectangularCuboidSeparatedVolumeMoments) {
  RectangularCuboid cube = unit_cell;
  Plane plane_up_y(Normal(0.0, 1.0, 0.0), 0.25);
  Plane plane_down_y(Normal(0.0, -1.0, 0.0), 0.25);
  PlanarSeparator recon_one_plane = PlanarSeparator::fromOnePlane(plane_up_y);

  // Cut cube by one plane reconstruction
  SeparatedMoments<VolumeMoments> from_cube_cut =
      getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          cube, recon_one_plane);
  double liquid_volume = 0.75;
  Pt liquid_centroid(0.0, -0.125, 0.0);
  double gas_volume = 0.25;
  Pt gas_centroid(0.0, 0.375, 0.0);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].volume(), liquid_volume);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].centroid()[0], liquid_centroid[0]);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].centroid()[1], liquid_centroid[1]);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].centroid()[2], liquid_centroid[2]);
  EXPECT_DOUBLE_EQ(from_cube_cut[1].volume(), gas_volume);
  EXPECT_DOUBLE_EQ(from_cube_cut[1].centroid()[0], gas_centroid[0]);
  EXPECT_DOUBLE_EQ(from_cube_cut[1].centroid()[1], gas_centroid[1]);
  EXPECT_DOUBLE_EQ(from_cube_cut[1].centroid()[2], gas_centroid[2]);

  // Cut cube by two plane reconstruction
  liquid_volume = 0.5;
  liquid_centroid = Pt(0.0, 0.0, 0.0);
  gas_volume = 0.5;
  gas_centroid = Pt(0.0, 0.0, 0.0);
  PlanarSeparator recon_two_plane =
      PlanarSeparator::fromTwoPlanes(plane_up_y, plane_down_y, 1.0);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, recon_two_plane);
  EXPECT_NEAR(from_cube_cut[0].volume(), liquid_volume, DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[2], liquid_centroid[2], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].volume(), gas_volume, DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], DBL_EPSILON);

  // Cut on shifted cube
  double shift_x = -4.8;
  double shift_y = 0.125;
  double shift_z = 3.28;
  cube.shift(shift_x, shift_y, shift_z);
  liquid_volume = 0.5;
  liquid_centroid = Pt(shift_x, 0.0, shift_z);
  gas_volume = 0.5;
  gas_centroid = Pt(shift_x, 2 * shift_y, shift_z);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, recon_two_plane);
  EXPECT_NEAR(from_cube_cut[0].volume(), liquid_volume, 1.0e-14);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], liquid_centroid[2], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].volume(), gas_volume, 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], 1.0e-14);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], 1.0e-14);

  cube.shift(-shift_x, -shift_y, -shift_z);
  Plane plane_up_x(Normal(1.0, 0.0, 0.0), 0.25);
  PlanarSeparator recon_one_plane_right =
      PlanarSeparator::fromOnePlane(plane_up_x);
  liquid_volume = 0.75;
  liquid_centroid = Pt(-0.125, 0.0, 0.0);
  gas_volume = 0.25;
  gas_centroid = Pt(0.375, 0.0, 0.0);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, recon_one_plane_right);
  EXPECT_NEAR(from_cube_cut[0].volume(), liquid_volume, DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[2], liquid_centroid[2], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].volume(), gas_volume, DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], DBL_EPSILON);

  Plane plane_45(Normal::normalized(1.0, 1.0, 0.0), 0.0);
  PlanarSeparator recon_45 = PlanarSeparator::fromOnePlane(plane_45);
  liquid_volume = 0.5;
  liquid_centroid = Pt(-0.5 / 3.0, -0.5 / 3.0, 0.0);
  gas_volume = 0.5;
  gas_centroid = Pt(0.5 / 3.0, 0.5 / 3.0, 0.0);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, recon_45);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].volume(), liquid_volume);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[2], liquid_centroid[2], DBL_EPSILON);
  EXPECT_DOUBLE_EQ(from_cube_cut[1].volume(), gas_volume);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], DBL_EPSILON);

  // Check extreme pure liquid/pure gas cases
  Plane one_plane_above_cell(Normal(0.0, 1.0, 0.0), 10.0);
  PlanarSeparator recon_one_above =
      PlanarSeparator::fromOnePlane(one_plane_above_cell);
  liquid_volume = 1.0;
  liquid_centroid = Pt(0.0, 0.0, 0.0);
  gas_volume = 0.0;
  gas_centroid = Pt(0.0, 0.0, 0.0);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, recon_one_above);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].volume(), liquid_volume);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[2], liquid_centroid[2], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].volume(), gas_volume, DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], DBL_EPSILON);

  one_plane_above_cell.normal() = Normal(0.0, -1.0, 0.0);
  one_plane_above_cell.distance() = (-10.0);
  recon_one_above[0] = one_plane_above_cell;
  liquid_volume = 0.0;
  liquid_centroid = Pt(0.0, 0.0, 0.0);
  gas_volume = 1.0;
  gas_centroid = Pt(0.0, 0.0, 0.0);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, recon_one_above);
  EXPECT_NEAR(from_cube_cut[0].volume(), liquid_volume, DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[2], liquid_centroid[2], DBL_EPSILON);
  EXPECT_DOUBLE_EQ(from_cube_cut[1].volume(), gas_volume);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], DBL_EPSILON);

  PlanarSeparator straddle_cell_liquid =
      PlanarSeparator::fromTwoPlanes(Plane(Normal(1.0, 0.0, 0.0), 2.0),
                                     Plane(Normal(-1.0, 0.0, 0.0), 2.0), 1.0);
  liquid_volume = 1.0;
  liquid_centroid = Pt(0.0, 0.0, 0.0);
  gas_volume = 0.0;
  gas_centroid = Pt(0.0, 0.0, 0.0);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, straddle_cell_liquid);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].volume(), liquid_volume);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[2], liquid_centroid[2], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].volume(), gas_volume, DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], DBL_EPSILON);

  PlanarSeparator straddle_cell_gas =
      PlanarSeparator::fromTwoPlanes(Plane(Normal(1.0, 0.0, 0.0), -2.0),
                                     Plane(Normal(-1.0, 0.0, 0.0), -2.0), 1.0);
  liquid_volume = 0.0;
  liquid_centroid = Pt(0.0, 0.0, 0.0);
  gas_volume = 1.0;
  gas_centroid = Pt(0.0, 0.0, 0.0);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, straddle_cell_gas);
  EXPECT_NEAR(from_cube_cut[0].volume(), liquid_volume, DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[2], liquid_centroid[2], DBL_EPSILON);
  EXPECT_DOUBLE_EQ(from_cube_cut[1].volume(), gas_volume);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], DBL_EPSILON);

  // Make sure flipping doesn't effect pure cells
  straddle_cell_liquid.flipCutting();
  liquid_volume = 1.0;
  liquid_centroid = Pt(0.0, 0.0, 0.0);
  gas_volume = 0.0;
  gas_centroid = Pt(0.0, 0.0, 0.0);
  from_cube_cut = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      cube, straddle_cell_liquid);
  EXPECT_DOUBLE_EQ(from_cube_cut[0].volume(), liquid_volume);
  EXPECT_NEAR(from_cube_cut[0].centroid()[0], liquid_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[1], liquid_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[0].centroid()[2], liquid_centroid[2], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].volume(), gas_volume, DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[0], gas_centroid[0], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[1], gas_centroid[1], DBL_EPSILON);
  EXPECT_NEAR(from_cube_cut[1].centroid()[2], gas_centroid[2], DBL_EPSILON);
}

TEST(GenericCutting, cutCuboidVolumeFraction) {
  Plane plane_up_y(Normal(0.0, 1.0, 0.0), 0.25);
  Plane plane_down_y(Normal(0.0, -1.0, 0.0), 0.25);
  PlanarSeparator recon_one_plane = PlanarSeparator::fromOnePlane(plane_up_y);

  // Cut cube by one plane reconstruction
  double cut_volume = getVolumeFraction(unit_cell, recon_one_plane);
  double liquid_volume = 0.75;
  EXPECT_DOUBLE_EQ(cut_volume, liquid_volume);

  // Cut cube by two plane reconstruction
  liquid_volume = 0.5;
  PlanarSeparator recon_two_plane =
      PlanarSeparator::fromTwoPlanes(plane_up_y, plane_down_y, 1.0);
  cut_volume = getVolumeFraction(unit_cell, recon_two_plane);
  EXPECT_DOUBLE_EQ(cut_volume, liquid_volume);

  Plane plane_up_x(Normal(1.0, 0.0, 0.0), 0.25);
  PlanarSeparator recon_one_plane_right =
      PlanarSeparator::fromOnePlane(plane_up_x);
  liquid_volume = 0.75;
  cut_volume = getVolumeFraction(unit_cell, recon_one_plane_right);
  EXPECT_DOUBLE_EQ(cut_volume, liquid_volume);

  Plane plane_45(Normal::normalized(1.0, 1.0, 0.0), 0.0);
  PlanarSeparator recon_45 = PlanarSeparator::fromOnePlane(plane_45);
  liquid_volume = 0.5;
  cut_volume = getVolumeFraction(unit_cell, recon_45);
  EXPECT_DOUBLE_EQ(cut_volume, liquid_volume);

  // Check extreme pure liquid/pure gas cases
  Plane one_plane_above_cell(Normal(0.0, 1.0, 0.0), 10.0);
  PlanarSeparator recon_one_above =
      PlanarSeparator::fromOnePlane(one_plane_above_cell);
  liquid_volume = 1.0;
  cut_volume = getVolumeFraction(unit_cell, recon_one_above);
  EXPECT_DOUBLE_EQ(cut_volume, liquid_volume);

  one_plane_above_cell.normal() = Normal(0.0, -1.0, 0.0);
  one_plane_above_cell.distance() = (-10.0);
  recon_one_above[0] = one_plane_above_cell;
  liquid_volume = 0.0;
  cut_volume = getVolumeFraction(unit_cell, recon_one_above);
  EXPECT_NEAR(cut_volume, liquid_volume, DBL_EPSILON);

  PlanarSeparator straddle_cell_liquid =
      PlanarSeparator::fromTwoPlanes(Plane(Normal(1.0, 0.0, 0.0), 2.0),
                                     Plane(Normal(-1.0, 0.0, 0.0), 2.0), 1.0);
  liquid_volume = 1.0;
  cut_volume = getVolumeFraction(unit_cell, straddle_cell_liquid);
  EXPECT_DOUBLE_EQ(cut_volume, liquid_volume);

  PlanarSeparator straddle_cell_gas =
      PlanarSeparator::fromTwoPlanes(Plane(Normal(1.0, 0.0, 0.0), -2.0),
                                     Plane(Normal(-1.0, 0.0, 0.0), -2.0), 1.0);
  liquid_volume = 0.0;
  cut_volume = getVolumeFraction(unit_cell, straddle_cell_gas);
  EXPECT_NEAR(cut_volume, liquid_volume, DBL_EPSILON);

  // Make sure flipping doesn't effect pure cells
  straddle_cell_liquid.flipCutting();
  liquid_volume = 1.0;
  cut_volume = getVolumeFraction(unit_cell, straddle_cell_liquid);
  EXPECT_DOUBLE_EQ(cut_volume, liquid_volume);
}

TEST(GenericCutting, cutTetByPlanarReconstructionForVolume) {
  Tet tet({Pt(0.0, 1.0, 0.0), Pt(0.0, 0.0, 1.0), Pt(1.0, 0.0, 0.0),
           Pt(0.0, 0.0, -1.0)});
  EXPECT_DOUBLE_EQ(tet.calculateVolume(), 0.5 * 2.0 / 3.0);

  Plane plane_up_y(Normal(0.0, 1.0, 0.0), 0.25);

  Pt midpoint_0 = Pt::fromEdgeIntersection(
      tet[0], plane_up_y.signedDistanceToPoint(tet[0]), tet[1],
      plane_up_y.signedDistanceToPoint(tet[1]));
  Pt midpoint_1 = Pt::fromEdgeIntersection(
      tet[0], plane_up_y.signedDistanceToPoint(tet[0]), tet[2],
      plane_up_y.signedDistanceToPoint(tet[2]));
  Pt midpoint_2 = Pt::fromEdgeIntersection(
      tet[0], plane_up_y.signedDistanceToPoint(tet[0]), tet[3],
      plane_up_y.signedDistanceToPoint(tet[3]));
  Tet upper_tet({tet[0], midpoint_0, midpoint_1, midpoint_2});

  // Case with plane_up_y
  PlanarSeparator recon_one_plane = PlanarSeparator::fromOnePlane(plane_up_y);
  // recon_one_plane.flipCutting();
  Volume cut_tet_volume =
      getNormalizedVolumeMoments<Volume>(tet, recon_one_plane);
  EXPECT_DOUBLE_EQ(tet.calculateVolume() - upper_tet.calculateVolume(),
                   cut_tet_volume);

  // Case where opposite side of plane is cut for
  recon_one_plane[0] = recon_one_plane[0].generateFlippedPlane();
  cut_tet_volume = getNormalizedVolumeMoments<Volume>(tet, recon_one_plane);
  EXPECT_DOUBLE_EQ(upper_tet.calculateVolume(), cut_tet_volume);

  // Test by using to cut cube by sheet
  // Cut cube by one plane reconstruction
  PlanarSeparator reconstruction_for_cube =
      PlanarSeparator::fromTwoPlanes(Plane(Normal(0.0, 1.0, 0.0), 0.25),
                                     Plane(Normal(0.0, -1.0, 0.0), 0.0), 1.0);
  VolumeMoments liquid_moments(0.0, Pt(0.0, 0.0, 0.0));
  double liquid_volume = 0.25;
  Volume cut_liquid_volume = 0.0;
  for (UnsignedIndex_t t = 0;
       t < unit_cell.getNumberOfSimplicesInDecomposition(); ++t) {
    cut_liquid_volume += getNormalizedVolumeMoments<Volume>(
        unit_cell.getSimplexFromDecomposition(t), reconstruction_for_cube);
  }
  EXPECT_DOUBLE_EQ(liquid_volume, cut_liquid_volume);
}

TEST(GenericCutting, localizedSeparatedVolumeMomentsTetByTet) {
  // Test when localizing lays entirely in surrounding
  Tet surrounding_tet({Pt(-2.0, 0.0, 0.0), Pt(0.0, 0.0, 2.0), Pt(2.0, 0.0, 0.0),
                       Pt(0.0, -2.0, 0.0)});
  Tet localizing_tet;
  for (UnsignedIndex_t v = 0; v < 4; ++v) {
    localizing_tet[v] = 0.5 * surrounding_tet[v];
  }
  PlanarSeparator far_above_plane =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 1.0, 0.0), 1000.0));
  PlanarLocalizer localizer = localizing_tet.getLocalizer();
  LocalizedSeparator local_sep(&localizer, &far_above_plane);
  SeparatedMoments<VolumeMoments> same_moments =
      getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          surrounding_tet, local_sep);
  Pt localizing_tet_centroid = localizing_tet.calculateCentroid();
  EXPECT_NEAR(same_moments[0].volume(), localizing_tet.calculateVolume(),
              1.0e-14);
  EXPECT_NEAR(same_moments[0].centroid()[0], localizing_tet_centroid[0],
              1.0e-14);
  EXPECT_NEAR(same_moments[0].centroid()[1], localizing_tet_centroid[1],
              1.0e-14);
  EXPECT_NEAR(same_moments[0].centroid()[2], localizing_tet_centroid[2],
              1.0e-14);

  PlanarSeparator far_below_plane =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 1.0, 0.0), -1000.0));
  local_sep = LocalizedSeparator(&localizer, &far_below_plane);
  same_moments = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      surrounding_tet, local_sep);
  EXPECT_NEAR(same_moments[1].volume(), localizing_tet.calculateVolume(),
              1.0e-14);
  EXPECT_NEAR(same_moments[1].centroid()[0], localizing_tet_centroid[0],
              1.0e-14);
  EXPECT_NEAR(same_moments[1].centroid()[1], localizing_tet_centroid[1],
              1.0e-14);
  EXPECT_NEAR(same_moments[1].centroid()[2], localizing_tet_centroid[2],
              1.0e-14);
}

TEST(GenericCutting, localizedSeparatedVolumeMomentsRectCuboidByRectCuboid) {
  RectangularCuboid corner_cube =
      RectangularCuboid::fromBoundingPts(Pt(0.0, 0.0, 0.0), Pt(0.5, 0.5, 0.5));
  PlanarSeparator far_above_plane =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 1.0, 0.0), 10000.0));
  PlanarLocalizer cc_localizer = corner_cube.getLocalizer();
  auto no_interface =
      getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          unit_cell, LocalizedSeparator(&cc_localizer, &far_above_plane));
  Pt corner_cube_centroid = corner_cube.calculateCentroid();
  EXPECT_NEAR(no_interface[0].volume(), corner_cube.calculateVolume(), 1.0e-14);
  EXPECT_NEAR(no_interface[0].centroid()[0], corner_cube_centroid[0], 1.0e-14);
  EXPECT_NEAR(no_interface[0].centroid()[1], corner_cube_centroid[1], 1.0e-14);
  EXPECT_NEAR(no_interface[0].centroid()[2], corner_cube_centroid[2], 1.0e-14);

  PlanarSeparator intersecting_interface = PlanarSeparator::fromTwoPlanes(
      Plane(Normal::normalized(-1.0, -1.0, -1.0), -0.1),
      Plane(Normal::normalized(1.0, 1.0, 1.0), -0.1), -1.0);

  auto correct_phase_moments =
      getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          corner_cube, intersecting_interface);
  auto localized_moments =
      getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          unit_cell,
          LocalizedSeparator(&cc_localizer, &intersecting_interface));
  EXPECT_NEAR(localized_moments[0].volume(), correct_phase_moments[0].volume(),
              1.0e-14);
  EXPECT_NEAR(localized_moments[0].centroid()[0],
              correct_phase_moments[0].centroid()[0], 1.0e-14);
  EXPECT_NEAR(localized_moments[0].centroid()[1],
              correct_phase_moments[0].centroid()[1], 1.0e-14);
  EXPECT_NEAR(localized_moments[0].centroid()[2],
              correct_phase_moments[0].centroid()[2], 1.0e-14);
  EXPECT_NEAR(localized_moments[1].volume(), correct_phase_moments[1].volume(),
              1.0e-14);
  EXPECT_NEAR(localized_moments[1].centroid()[0],
              correct_phase_moments[1].centroid()[0], 1.0e-14);
  EXPECT_NEAR(localized_moments[1].centroid()[1],
              correct_phase_moments[1].centroid()[1], 1.0e-14);
  EXPECT_NEAR(localized_moments[1].centroid()[2],
              correct_phase_moments[1].centroid()[2], 1.0e-14);
}

TEST(GenericCutting, extremaLocalizerCuts) {
  RectangularCuboid small_cell = RectangularCuboid::fromBoundingPts(
      Pt(-0.25, -0.25, -0.25), Pt(0.25, 0.25, 0.25));
  RectangularCuboid bigger_cell = unit_cell;
  PlanarSeparator flat =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 1.0, 0.0), 0.0));
  PlanarLocalizer bigcell_localizer = bigger_cell.getLocalizer();
  auto phase_moments =
      getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          small_cell, LocalizedSeparator(&bigcell_localizer, &flat));
  auto correct_phase_moments =
      getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(small_cell,
                                                                  flat);
  EXPECT_NEAR(phase_moments[0].volume(), correct_phase_moments[0].volume(),
              1.0e-14);
  EXPECT_NEAR(phase_moments[0].centroid()[0],
              correct_phase_moments[0].centroid()[0], 1.0e-14);
  EXPECT_NEAR(phase_moments[0].centroid()[1],
              correct_phase_moments[0].centroid()[1], 1.0e-14);
  EXPECT_NEAR(phase_moments[0].centroid()[2],
              correct_phase_moments[0].centroid()[2], 1.0e-14);
  EXPECT_NEAR(phase_moments[1].volume(), correct_phase_moments[1].volume(),
              1.0e-14);
  EXPECT_NEAR(phase_moments[1].centroid()[0],
              correct_phase_moments[1].centroid()[0], 1.0e-14);
  EXPECT_NEAR(phase_moments[1].centroid()[1],
              correct_phase_moments[1].centroid()[1], 1.0e-14);
  EXPECT_NEAR(phase_moments[1].centroid()[2],
              correct_phase_moments[1].centroid()[2], 1.0e-14);

  bigger_cell.shift(-10.0, 0.0, 0.0);
  bigcell_localizer = bigger_cell.getLocalizer();
  phase_moments = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      small_cell, LocalizedSeparator(&bigcell_localizer, &flat));
  EXPECT_NEAR(phase_moments[0].volume(), 0.0, 1.0e-14);
  EXPECT_NEAR(phase_moments[0].centroid()[0], 0.0, 1.0e-14);
  EXPECT_NEAR(phase_moments[0].centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(phase_moments[0].centroid()[2], 0.0, 1.0e-14);
  EXPECT_NEAR(phase_moments[1].volume(), 0.0, 1.0e-14);
  EXPECT_NEAR(phase_moments[1].centroid()[0], 0.0, 1.0e-14);
  EXPECT_NEAR(phase_moments[1].centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(phase_moments[1].centroid()[2], 0.0, 1.0e-14);
}

TEST(GenericCutting, simpleLinkage) {
  RectangularCuboid super_cell = RectangularCuboid::fromBoundingPts(
      Pt(-5.5, -0.5, -0.5), Pt(-1.5, 0.5, 0.5));
  RectangularCuboid cell[4];
  PlanarLocalizer localizers[4];
  PlanarSeparator separators[4];
  LocalizedSeparatorLink link_combined[4];

  cell[0] = RectangularCuboid::fromBoundingPts(Pt(-5.5, -0.5, -0.5),
                                               Pt(-4.5, 0.5, 0.5));
  cell[1] = RectangularCuboid::fromBoundingPts(Pt(-4.5, -0.5, -0.5),
                                               Pt(-3.5, 0.5, 0.5));
  cell[2] = RectangularCuboid::fromBoundingPts(Pt(-3.5, -0.5, -0.5),
                                               Pt(-2.5, 0.5, 0.5));
  cell[3] = RectangularCuboid::fromBoundingPts(Pt(-2.5, -0.5, -0.5),
                                               Pt(-1.5, 0.5, 0.5));

  localizers[0] = cell[0].getLocalizer();
  localizers[1] = cell[1].getLocalizer();
  localizers[2] = cell[2].getLocalizer();
  localizers[3] = cell[3].getLocalizer();

  separators[0] =
      PlanarSeparator::fromOnePlane(Plane(Normal(-1.0, 0.0, 0.0), 5.0));
  separators[1] =
      PlanarSeparator::fromOnePlane(Plane(Normal(-1.0, 0.0, 0.0), 5.0));
  separators[2] =
      PlanarSeparator::fromOnePlane(Plane(Normal(1.0, 0.0, 0.0), -2.0));
  separators[3] =
      PlanarSeparator::fromOnePlane(Plane(Normal(1.0, 0.0, 0.0), -2.0));
  for (UnsignedIndex_t n = 0; n < 4; ++n) {
    link_combined[n] = LocalizedSeparatorLink(&localizers[n], &separators[n]);
    link_combined[n].setId(n);
    if (n != 0) {
      link_combined[n].setEdgeConnectivity(0, &link_combined[n - 1]);
    }
    if (n != 3) {
      link_combined[n].setEdgeConnectivity(1, &link_combined[n + 1]);
    }
  }

  auto summed_moments =
      getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          super_cell, link_combined[2]);

  EXPECT_NEAR(summed_moments[0].volume(), 3.0, 1.0e-14);
  EXPECT_NEAR(summed_moments[1].volume(), 1.0, 1.0e-14);
  EXPECT_NEAR(summed_moments[0].centroid()[0], -3.5, 1.0e-14);
  EXPECT_NEAR(summed_moments[0].centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(summed_moments[0].centroid()[2], 0.0, 1.0e-14);
  EXPECT_NEAR(summed_moments[1].centroid()[0], -3.5, 1.0e-14);
  EXPECT_NEAR(summed_moments[1].centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(summed_moments[1].centroid()[2], 0.0, 1.0e-14);
}

TEST(GenericCutting, getVolumeMomentsLocalizedSeparatorLink) {
  RectangularCuboid super_cell = RectangularCuboid::fromBoundingPts(
      Pt(-1.5, -1.5, -0.5), Pt(1.5, 1.5, 0.5));
  RectangularCuboid cell[3][3];
  PlanarLocalizer localizers[3][3];
  PlanarSeparator separators[3][3];
  LocalizedSeparatorLink link_combined[3][3];
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      cell[j][i] = unit_cell;
      cell[j][i].shift(static_cast<double>(i - 1), static_cast<double>(j - 1),
                       0.0);
      localizers[j][i] = cell[j][i].getLocalizer();
    }
  }
  Normal plane_normals[3][3];
  plane_normals[0][0] = Normal(-1.0, -1.0, 0.0);
  plane_normals[0][2] = Normal(1.0, -1.0, 0.0);
  plane_normals[2][0] = Normal(-1.0, 1.0, 0.0);
  plane_normals[2][2] = Normal(1.0, 1.0, 0.0);
  plane_normals[0][1] = Normal(0.0, -1.0, 0.0);
  plane_normals[2][1] = Normal(0.0, 1.0, 0.0);
  plane_normals[1][0] = Normal(-1.0, 0.0, 0.0);
  plane_normals[1][2] = Normal(1.0, 0.0, 0.0);
  plane_normals[1][1] = Normal(0.0, 0.0, 0.0);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      plane_normals[j][i].normalize();
      separators[j][i] = PlanarSeparator::fromOnePlane(
          Plane(plane_normals[j][i],
                plane_normals[j][i] * (cell[j][i].calculateCentroid())));
      link_combined[j][i] =
          LocalizedSeparatorLink(&localizers[j][i], &separators[j][i]);
    }
  }
  setToPurePhaseReconstruction(1.0, &separators[1][1]);

  LocalizedSeparatorLink* ptr;
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      link_combined[j][i].setId(static_cast<UnsignedIndex_t>(j * 3 + i));
      ptr = i - 1 >= 0 ? &link_combined[j][i - 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(0, ptr);
      ptr = i + 1 < 3 ? &link_combined[j][i + 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(1, ptr);
      ptr = j - 1 >= 0 ? &link_combined[j - 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(2, ptr);
      ptr = j + 1 < 3 ? &link_combined[j + 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(3, ptr);
    }
  }

  SeparatedMoments<VolumeMoments> correct_volume_moments;
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      correct_volume_moments +=
          getVolumeMoments<SeparatedMoments<VolumeMoments>>(cell[j][i],
                                                            separators[j][i]);
    }
  }
  EXPECT_NEAR(correct_volume_moments[0].volume(), 5.0, 1.0e-14);
  EXPECT_NEAR(correct_volume_moments[1].volume(), 4.0, 1.0e-14);
  EXPECT_NEAR(correct_volume_moments[0].centroid()[0], 0.0, 1.0e-14);
  EXPECT_NEAR(correct_volume_moments[0].centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(correct_volume_moments[0].centroid()[2], 0.0, 1.0e-14);
  EXPECT_NEAR(correct_volume_moments[1].centroid()[0], 0.0, 1.0e-14);
  EXPECT_NEAR(correct_volume_moments[1].centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(correct_volume_moments[1].centroid()[2], 0.0, 1.0e-14);

  auto volume_moments = getVolumeMoments<SeparatedMoments<VolumeMoments>>(
      super_cell, link_combined[0][1]);
  EXPECT_NEAR(volume_moments[0].volume(), 5.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[1].volume(), 4.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[0].centroid()[0], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[0].centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[0].centroid()[2], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[1].centroid()[0], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[1].centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[1].centroid()[2], 0.0, 1.0e-14);

  // Should be agnostic to starting location
  volume_moments = getVolumeMoments<SeparatedMoments<VolumeMoments>>(
      super_cell, link_combined[2][2]);
  EXPECT_NEAR(volume_moments[0].volume(), 5.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[1].volume(), 4.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[0].centroid()[0], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[0].centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[0].centroid()[2], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[1].centroid()[0], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[1].centroid()[1], 0.0, 1.0e-14);
  EXPECT_NEAR(volume_moments[1].centroid()[2], 0.0, 1.0e-14);
}

TEST(GenericCutting, getVolumeMomentsListed) {
  RectangularCuboid super_cell = RectangularCuboid::fromBoundingPts(
      Pt(-1.5, -1.5, -0.5), Pt(1.5, 1.5, 0.5));
  RectangularCuboid cell[3][3];
  PlanarLocalizer localizers[3][3];
  PlanarSeparator separators[3][3];
  LocalizedSeparatorLink link_combined[3][3];
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      cell[j][i] = unit_cell;
      cell[j][i].shift(static_cast<double>(i - 1), static_cast<double>(j - 1),
                       0.0);
      localizers[j][i] = cell[j][i].getLocalizer();
    }
  }
  Normal plane_normals[3][3];
  plane_normals[0][0] = Normal(-1.0, -1.0, 0.0);
  plane_normals[0][2] = Normal(1.0, -1.0, 0.0);
  plane_normals[2][0] = Normal(-1.0, 1.0, 0.0);
  plane_normals[2][2] = Normal(1.0, 1.0, 0.0);
  plane_normals[0][1] = Normal(0.0, -1.0, 0.0);
  plane_normals[2][1] = Normal(0.0, 1.0, 0.0);
  plane_normals[1][0] = Normal(-1.0, 0.0, 0.0);
  plane_normals[1][2] = Normal(1.0, 0.0, 0.0);
  plane_normals[1][1] = Normal(0.0, 0.0, 0.0);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      plane_normals[j][i].normalize();
      separators[j][i] = PlanarSeparator::fromOnePlane(
          Plane(plane_normals[j][i],
                plane_normals[j][i] * (cell[j][i].calculateCentroid())));
      link_combined[j][i] =
          LocalizedSeparatorLink(&localizers[j][i], &separators[j][i]);
    }
  }

  setToPurePhaseReconstruction(1.0, &separators[1][1]);

  LocalizedSeparatorLink* ptr;
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      link_combined[j][i].setId(static_cast<UnsignedIndex_t>(j * 3 + i));
      ptr = i - 1 >= 0 ? &link_combined[j][i - 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(0, ptr);
      ptr = i + 1 < 3 ? &link_combined[j][i + 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(1, ptr);
      ptr = j - 1 >= 0 ? &link_combined[j - 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(2, ptr);
      ptr = j + 1 < 3 ? &link_combined[j + 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(3, ptr);
    }
  }

  SeparatedMoments<VolumeMoments> correct_volume_moments[3][3];
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      correct_volume_moments[j][i] =
          getVolumeMoments<SeparatedMoments<VolumeMoments>>(cell[j][i],
                                                            separators[j][i]);
    }
  }
  auto volume_moments = getVolumeMoments<
      AccumulatedVolumeMoments<SeparatedMoments<VolumeMoments>>>(
      super_cell, link_combined[0][1]);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      auto ind = static_cast<UnsignedIndex_t>(j * 3 + i);
      EXPECT_NEAR(volume_moments[ind][0].volume(),
                  correct_volume_moments[j][i][0].volume(), 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][1].volume(),
                  correct_volume_moments[j][i][1].volume(), 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[0],
                  correct_volume_moments[j][i][0].centroid()[0], 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[1],
                  correct_volume_moments[j][i][0].centroid()[1], 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[2],
                  correct_volume_moments[j][i][0].centroid()[2], 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[0],
                  correct_volume_moments[j][i][1].centroid()[0], 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[1],
                  correct_volume_moments[j][i][1].centroid()[1], 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[2],
                  correct_volume_moments[j][i][1].centroid()[2], 1.0e-15);
    }
  }
  // Should be agnostic to starting location.
  volume_moments = getVolumeMoments<
      AccumulatedVolumeMoments<SeparatedMoments<VolumeMoments>>>(
      super_cell, link_combined[2][2]);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      auto ind = static_cast<UnsignedIndex_t>(j * 3 + i);
      EXPECT_NEAR(volume_moments[ind][0].volume(),
                  correct_volume_moments[j][i][0].volume(), 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][1].volume(),
                  correct_volume_moments[j][i][1].volume(), 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[0],
                  correct_volume_moments[j][i][0].centroid()[0], 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[1],
                  correct_volume_moments[j][i][0].centroid()[1], 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[2],
                  correct_volume_moments[j][i][0].centroid()[2], 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[0],
                  correct_volume_moments[j][i][1].centroid()[0], 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[1],
                  correct_volume_moments[j][i][1].centroid()[1], 1.0e-15);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[2],
                  correct_volume_moments[j][i][1].centroid()[2], 1.0e-15);
    }
  }
}

TEST(GenericCutting, getVolumeMomentsNonCubicLocalizer) {
  RectangularCuboid super_cell = RectangularCuboid::fromBoundingPts(
      Pt(-1.5, -1.5, -0.5), Pt(1.5, 1.5, 0.5));
  RectangularCuboid cell[3][3];
  PlanarLocalizer localizers[3][3];
  PlanarSeparator separators[3][3];
  LocalizedSeparatorLink link_combined[3][3];
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      cell[j][i] = unit_cell;
      cell[j][i].shift(static_cast<double>(i - 1), static_cast<double>(j - 1),
                       0.0);
      localizers[j][i] = cell[j][i].getLocalizer();
    }
  }
  Normal plane_normals[3][3];
  plane_normals[0][0] = Normal(-1.0, -1.0, 0.0);
  plane_normals[0][2] = Normal(1.0, -1.0, 0.0);
  plane_normals[2][0] = Normal(-1.0, 1.0, 0.0);
  plane_normals[2][2] = Normal(1.0, 1.0, 0.0);
  plane_normals[0][1] = Normal(0.0, -1.0, 0.0);
  plane_normals[2][1] = Normal(0.0, 1.0, 0.0);
  plane_normals[1][0] = Normal(-1.0, 0.0, 0.0);
  plane_normals[1][2] = Normal(1.0, 0.0, 0.0);
  plane_normals[1][1] = Normal(0.0, 0.0, 0.0);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      plane_normals[j][i].normalize();
      plane_normals[j][i] = -plane_normals[j][i];
      localizers[j][i].addPlane(
          Plane(plane_normals[j][i],
                plane_normals[j][i] * (cell[j][i].calculateCentroid())));
      separators[j][i] =
          PlanarSeparator::fromOnePlane(Plane(Normal(-1.0, 0.0, 0.0), 0.0));
      link_combined[j][i] =
          LocalizedSeparatorLink(&localizers[j][i], &separators[j][i]);
    }
  }
  (localizers[1][1])[6].distance() = (-1000000.00);

  LocalizedSeparatorLink* ptr;
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      link_combined[j][i].setId(static_cast<UnsignedIndex_t>(j * 3 + i));
      ptr = i - 1 >= 0 ? &link_combined[j][i - 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(0, ptr);
      ptr = i + 1 < 3 ? &link_combined[j][i + 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(1, ptr);
      ptr = j - 1 >= 0 ? &link_combined[j - 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(2, ptr);
      ptr = j + 1 < 3 ? &link_combined[j + 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(3, ptr);
    }
  }

  SeparatedMoments<VolumeMoments> correct_volume_moments[3][3];
  double ss = 7.0 / 6.0;
  correct_volume_moments[0][0] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.00, Pt(-1000.0, -1000.0, -1000.0)),
      VolumeMoments(0.50, Pt(-ss, -ss, 0.0)));
  correct_volume_moments[0][1] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.25, Pt(0.25, -1.25, 0.0)),
      VolumeMoments(0.25, Pt(-0.25, -1.25, 0.0)));
  correct_volume_moments[0][2] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.50, Pt(ss, -ss, 0.0)),
      VolumeMoments(0.00, Pt(-1000.0, -1000.0, -1000.0)));
  correct_volume_moments[1][0] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.00, Pt(-1000.0, -1000.0, -1000.0)),
      VolumeMoments(0.50, Pt(-1.25, 0.0, 0.0)));
  correct_volume_moments[1][1] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.00, Pt(-1000.0, -1000.0, -1000.0)),
      VolumeMoments(0.00, Pt(-1000.0, -1000.0, -1000.0)));
  correct_volume_moments[1][2] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.50, Pt(1.25, 0.0, 0.0)),
      VolumeMoments(0.00, Pt(-1000.0, -1000.0, -1000.0)));
  correct_volume_moments[2][0] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.00, Pt(-1000.0, -1000.0, -1000.0)),
      VolumeMoments(0.50, Pt(-ss, ss, 0.0)));
  correct_volume_moments[2][1] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.25, Pt(0.25, 1.25, 0.0)),
      VolumeMoments(0.25, Pt(-0.25, 1.25, 0.0)));
  correct_volume_moments[2][2] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.50, Pt(ss, ss, 0.0)),
      VolumeMoments(0.00, Pt(-1000.0, -1000.0, -1000.0)));

  SeparatedMoments<VolumeMoments> whole_moments;
  double volume = {0.0};
  Pt centroid(0.0, 0.0, 0.0);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      correct_volume_moments[j][i].multiplyByVolume();
      whole_moments += correct_volume_moments[j][i];
    }
  }
  whole_moments.normalizeByVolume();
  // Should be agnostic to starting location
  auto volume_moments =
      getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          super_cell, link_combined[2][2]);
  EXPECT_NEAR(volume_moments[0].volume(), whole_moments[0].volume(), 1.0e-14);
  EXPECT_NEAR(volume_moments[1].volume(), whole_moments[1].volume(), 1.0e-14);
  EXPECT_NEAR(volume_moments[0].centroid()[0], whole_moments[0].centroid()[0],
              1.0e-14);
  EXPECT_NEAR(volume_moments[0].centroid()[1], whole_moments[0].centroid()[1],
              1.0e-14);
  EXPECT_NEAR(volume_moments[0].centroid()[2], whole_moments[0].centroid()[2],
              1.0e-14);
  EXPECT_NEAR(volume_moments[1].centroid()[0], whole_moments[1].centroid()[0],
              1.0e-14);
  EXPECT_NEAR(volume_moments[1].centroid()[1], whole_moments[1].centroid()[1],
              1.0e-14);
  EXPECT_NEAR(volume_moments[1].centroid()[2], whole_moments[1].centroid()[2],
              1.0e-14);

  // Should be agnostic to starting location
  volume_moments = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
      super_cell, link_combined[1][0]);
  EXPECT_NEAR(volume_moments[0].volume(), whole_moments[0].volume(), 1.0e-14);
  EXPECT_NEAR(volume_moments[1].volume(), whole_moments[1].volume(), 1.0e-14);
  EXPECT_NEAR(volume_moments[0].centroid()[0], whole_moments[0].centroid()[0],
              1.0e-14);
  EXPECT_NEAR(volume_moments[0].centroid()[1], whole_moments[0].centroid()[1],
              1.0e-14);
  EXPECT_NEAR(volume_moments[0].centroid()[2], whole_moments[0].centroid()[2],
              1.0e-14);
  EXPECT_NEAR(volume_moments[1].centroid()[0], whole_moments[1].centroid()[0],
              1.0e-14);
  EXPECT_NEAR(volume_moments[1].centroid()[1], whole_moments[1].centroid()[1],
              1.0e-14);
  EXPECT_NEAR(volume_moments[1].centroid()[2], whole_moments[1].centroid()[2],
              1.0e-14);
}

TEST(GenericCutting, getVolumeMomentsNonCubicLocalizerAccumulated) {
  RectangularCuboid super_cell = RectangularCuboid::fromBoundingPts(
      Pt(-1.5, -1.5, -0.5), Pt(1.5, 1.5, 0.5));
  RectangularCuboid cell[3][3];
  PlanarLocalizer localizers[3][3];
  PlanarSeparator separators[3][3];
  LocalizedSeparatorLink link_combined[3][3];
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      cell[j][i] = unit_cell;
      cell[j][i].shift(static_cast<double>(i - 1), static_cast<double>(j - 1),
                       0.0);
      localizers[j][i] = cell[j][i].getLocalizer();
    }
  }
  Normal plane_normals[3][3];
  plane_normals[0][0] = Normal(-1.0, -1.0, 0.0);
  plane_normals[0][2] = Normal(1.0, -1.0, 0.0);
  plane_normals[2][0] = Normal(-1.0, 1.0, 0.0);
  plane_normals[2][2] = Normal(1.0, 1.0, 0.0);
  plane_normals[0][1] = Normal(0.0, -1.0, 0.0);
  plane_normals[2][1] = Normal(0.0, 1.0, 0.0);
  plane_normals[1][0] = Normal(-1.0, 0.0, 0.0);
  plane_normals[1][2] = Normal(1.0, 0.0, 0.0);
  plane_normals[1][1] = Normal(0.0, 0.0, 0.0);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      plane_normals[j][i].normalize();
      plane_normals[j][i] = -plane_normals[j][i];
      localizers[j][i].addPlane(
          Plane(plane_normals[j][i],
                plane_normals[j][i] * (cell[j][i].calculateCentroid())));
      separators[j][i] =
          PlanarSeparator::fromOnePlane(Plane(Normal(-1.0, 0.0, 0.0), 0.0));
      link_combined[j][i] =
          LocalizedSeparatorLink(&localizers[j][i], &separators[j][i]);
    }
  }
  (localizers[1][1])[6].distance() = (-1000000.00);

  SeparatedMoments<VolumeMoments> correct_volume_moments[3][3];
  double ss = 7.0 / 6.0;
  correct_volume_moments[0][0] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.50, Pt(-ss, -ss, 0.0)));
  correct_volume_moments[0][1] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.25, Pt(0.25, -1.25, 0.0)),
      VolumeMoments(0.25, Pt(-0.25, -1.25, 0.0)));
  correct_volume_moments[0][2] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.50, Pt(ss, -ss, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));
  correct_volume_moments[1][0] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.50, Pt(-1.25, 0.0, 0.0)));
  correct_volume_moments[1][1] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));
  correct_volume_moments[1][2] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.50, Pt(1.25, 0.0, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));
  correct_volume_moments[2][0] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.50, Pt(-ss, ss, 0.0)));
  correct_volume_moments[2][1] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.25, Pt(0.25, 1.25, 0.0)),
      VolumeMoments(0.25, Pt(-0.25, 1.25, 0.0)));
  correct_volume_moments[2][2] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.50, Pt(ss, ss, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));

  // Check individual cutting by these LocalizedSeparators
  SeparatedMoments<VolumeMoments> individual_volume_moments[9];
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      individual_volume_moments[j * 3 + i] =
          getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
              super_cell, link_combined[j][i]);
    }
  }
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      auto ind = static_cast<UnsignedIndex_t>(j * 3 + i);
      EXPECT_NEAR(individual_volume_moments[ind][0].volume(),
                  correct_volume_moments[j][i][0].volume(), 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].volume(),
                  correct_volume_moments[j][i][1].volume(), 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][0].centroid()[0],
                  correct_volume_moments[j][i][0].centroid()[0], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][0].centroid()[1],
                  correct_volume_moments[j][i][0].centroid()[1], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][0].centroid()[2],
                  correct_volume_moments[j][i][0].centroid()[2], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].centroid()[0],
                  correct_volume_moments[j][i][1].centroid()[0], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].centroid()[1],
                  correct_volume_moments[j][i][1].centroid()[1], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].centroid()[2],
                  correct_volume_moments[j][i][1].centroid()[2], 1.0e-14);
    }
  }

  LocalizedSeparatorLink* ptr;
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      link_combined[j][i].setId(static_cast<UnsignedIndex_t>(j * 3 + i));
      ptr = i - 1 >= 0 ? &link_combined[j][i - 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(0, ptr);
      ptr = i + 1 < 3 ? &link_combined[j][i + 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(1, ptr);
      ptr = j - 1 >= 0 ? &link_combined[j - 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(2, ptr);
      ptr = j + 1 < 3 ? &link_combined[j + 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(3, ptr);
    }
  }

  // Should be agnostic to starting location
  auto volume_moments = getNormalizedVolumeMoments<
      AccumulatedVolumeMoments<SeparatedMoments<VolumeMoments>>>(
      super_cell, link_combined[0][0]);
  EXPECT_EQ(volume_moments.size(), 9);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      auto ind = static_cast<UnsignedIndex_t>(j * 3 + i);
      EXPECT_NEAR(volume_moments[ind][0].volume(),
                  correct_volume_moments[j][i][0].volume(), 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].volume(),
                  correct_volume_moments[j][i][1].volume(), 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[0],
                  correct_volume_moments[j][i][0].centroid()[0], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[1],
                  correct_volume_moments[j][i][0].centroid()[1], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[2],
                  correct_volume_moments[j][i][0].centroid()[2], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[0],
                  correct_volume_moments[j][i][1].centroid()[0], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[1],
                  correct_volume_moments[j][i][1].centroid()[1], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[2],
                  correct_volume_moments[j][i][1].centroid()[2], 1.0e-14);
    }
  }

  // Should be agnostic to starting location
  volume_moments = getNormalizedVolumeMoments<
      AccumulatedVolumeMoments<SeparatedMoments<VolumeMoments>>>(
      super_cell, link_combined[2][2]);
  EXPECT_EQ(volume_moments.size(), 9);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      auto ind = static_cast<UnsignedIndex_t>(j * 3 + i);
      EXPECT_NEAR(volume_moments[ind][0].volume(),
                  correct_volume_moments[j][i][0].volume(), 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].volume(),
                  correct_volume_moments[j][i][1].volume(), 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[0],
                  correct_volume_moments[j][i][0].centroid()[0], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[1],
                  correct_volume_moments[j][i][0].centroid()[1], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[2],
                  correct_volume_moments[j][i][0].centroid()[2], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[0],
                  correct_volume_moments[j][i][1].centroid()[0], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[1],
                  correct_volume_moments[j][i][1].centroid()[1], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[2],
                  correct_volume_moments[j][i][1].centroid()[2], 1.0e-14);
    }
  }
}

TEST(GenericCutting, getVolumeMomentsNonCubicLocalizerTagged) {
  RectangularCuboid super_cell = RectangularCuboid::fromBoundingPts(
      Pt(-1.5, -1.5, -0.5), Pt(1.5, 1.5, 0.5));
  RectangularCuboid cell[3][3];
  PlanarLocalizer localizers[3][3];
  PlanarSeparator separators[3][3];
  LocalizedSeparatorLink link_combined[3][3];
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      cell[j][i] = unit_cell;
      cell[j][i].shift(static_cast<double>(i - 1), static_cast<double>(j - 1),
                       0.0);
      localizers[j][i] = cell[j][i].getLocalizer();
    }
  }
  Normal plane_normals[3][3];
  plane_normals[0][0] = Normal(-1.0, -1.0, 0.0);
  plane_normals[0][2] = Normal(1.0, -1.0, 0.0);
  plane_normals[2][0] = Normal(-1.0, 1.0, 0.0);
  plane_normals[2][2] = Normal(1.0, 1.0, 0.0);
  plane_normals[0][1] = Normal(0.0, -1.0, 0.0);
  plane_normals[2][1] = Normal(0.0, 1.0, 0.0);
  plane_normals[1][0] = Normal(-1.0, 0.0, 0.0);
  plane_normals[1][2] = Normal(1.0, 0.0, 0.0);
  plane_normals[1][1] = Normal(0.0, 0.0, 0.0);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      plane_normals[j][i].normalize();
      plane_normals[j][i] = -plane_normals[j][i];
      localizers[j][i].addPlane(
          Plane(plane_normals[j][i],
                plane_normals[j][i] * (cell[j][i].calculateCentroid())));
      separators[j][i] =
          PlanarSeparator::fromOnePlane(Plane(Normal(-1.0, 0.0, 0.0), 0.0));
      link_combined[j][i] =
          LocalizedSeparatorLink(&localizers[j][i], &separators[j][i]);
    }
  }
  (localizers[1][1])[6].distance() = (-1000000.00);

  SeparatedMoments<VolumeMoments> correct_volume_moments[3][3];
  double ss = 7.0 / 6.0;
  correct_volume_moments[0][0] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.50, Pt(-ss, -ss, 0.0)));
  correct_volume_moments[0][1] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.25, Pt(0.25, -1.25, 0.0)),
      VolumeMoments(0.25, Pt(-0.25, -1.25, 0.0)));
  correct_volume_moments[0][2] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.50, Pt(ss, -ss, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));
  correct_volume_moments[1][0] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.50, Pt(-1.25, 0.0, 0.0)));
  correct_volume_moments[1][1] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));
  correct_volume_moments[1][2] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.50, Pt(1.25, 0.0, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));
  correct_volume_moments[2][0] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.50, Pt(-ss, ss, 0.0)));
  correct_volume_moments[2][1] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.25, Pt(0.25, 1.25, 0.0)),
      VolumeMoments(0.25, Pt(-0.25, 1.25, 0.0)));
  correct_volume_moments[2][2] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.50, Pt(ss, ss, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));

  // Check individual cutting by these LocalizedSeparators
  SeparatedMoments<VolumeMoments> individual_volume_moments[9];
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      individual_volume_moments[j * 3 + i] =
          getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
              super_cell, link_combined[j][i]);
    }
  }
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      auto ind = static_cast<UnsignedIndex_t>(j * 3 + i);
      EXPECT_NEAR(individual_volume_moments[ind][0].volume(),
                  correct_volume_moments[j][i][0].volume(), 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].volume(),
                  correct_volume_moments[j][i][1].volume(), 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][0].centroid()[0],
                  correct_volume_moments[j][i][0].centroid()[0], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][0].centroid()[1],
                  correct_volume_moments[j][i][0].centroid()[1], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][0].centroid()[2],
                  correct_volume_moments[j][i][0].centroid()[2], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].centroid()[0],
                  correct_volume_moments[j][i][1].centroid()[0], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].centroid()[1],
                  correct_volume_moments[j][i][1].centroid()[1], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].centroid()[2],
                  correct_volume_moments[j][i][1].centroid()[2], 1.0e-14);
    }
  }

  LocalizedSeparatorLink* ptr;
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      link_combined[j][i].setId(
          static_cast<UnsignedIndex_t>(8 * j * 3 + i));  // Non-contigous tag
      ptr = i - 1 >= 0 ? &link_combined[j][i - 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(0, ptr);
      ptr = i + 1 < 3 ? &link_combined[j][i + 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(1, ptr);
      ptr = j - 1 >= 0 ? &link_combined[j - 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(2, ptr);
      ptr = j + 1 < 3 ? &link_combined[j + 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(3, ptr);
    }
  }

  // Should be agnostic to starting location
  auto volume_moments = getNormalizedVolumeMoments<
      TaggedAccumulatedVolumeMoments<SeparatedMoments<VolumeMoments>>>(
      super_cell, link_combined[0][0]);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      // Lookup via known tag value
      auto tag = static_cast<UnsignedIndex_t>(8 * j * 3 + i);
      EXPECT_NEAR(volume_moments[tag][0].volume(),
                  correct_volume_moments[j][i][0].volume(), 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][1].volume(),
                  correct_volume_moments[j][i][1].volume(), 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][0].centroid()[0],
                  correct_volume_moments[j][i][0].centroid()[0], 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][0].centroid()[1],
                  correct_volume_moments[j][i][0].centroid()[1], 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][0].centroid()[2],
                  correct_volume_moments[j][i][0].centroid()[2], 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][1].centroid()[0],
                  correct_volume_moments[j][i][1].centroid()[0], 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][1].centroid()[1],
                  correct_volume_moments[j][i][1].centroid()[1], 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][1].centroid()[2],
                  correct_volume_moments[j][i][1].centroid()[2], 1.0e-14);
    }
  }

  // Should be agnostic to starting location
  volume_moments = getNormalizedVolumeMoments<
      TaggedAccumulatedVolumeMoments<SeparatedMoments<VolumeMoments>>>(
      super_cell, link_combined[2][2]);
  for (UnsignedIndex_t elem = 0; elem < volume_moments.size(); ++elem) {
    auto tag = volume_moments.getTagForIndex(elem);
    int j = static_cast<int>(tag / 24);
    int i = static_cast<int>(tag) - j * 24;
    // Get tag value for entry and compare to correct moments
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[0].volume(),
                correct_volume_moments[j][i][0].volume(), 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[1].volume(),
                correct_volume_moments[j][i][1].volume(), 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[0].centroid()[0],
                correct_volume_moments[j][i][0].centroid()[0], 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[0].centroid()[1],
                correct_volume_moments[j][i][0].centroid()[1], 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[0].centroid()[2],
                correct_volume_moments[j][i][0].centroid()[2], 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[1].centroid()[0],
                correct_volume_moments[j][i][1].centroid()[0], 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[1].centroid()[1],
                correct_volume_moments[j][i][1].centroid()[1], 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[1].centroid()[2],
                correct_volume_moments[j][i][1].centroid()[2], 1.0e-14);
  }
}

TEST(GenericCutting, getVolumeMomentsNonCubicLocalizerListedTri) {
  Polygon poly_super_cell;
  poly_super_cell.addVertex(Pt(1.5, -1.5, 0.0));
  poly_super_cell.addVertex(Pt(1.5, 1.5, 0.0));
  poly_super_cell.addVertex(Pt(-1.5, 1.5, 0.0));
  poly_super_cell.addVertex(Pt(-1.5, -1.5, 0.0));
  poly_super_cell.calculateAndSetPlaneOfExistence();
  DividedPolygon super_cell = DividedPolygon::fromPolygon(poly_super_cell);
  RectangularCuboid cell[3][3];
  PlanarLocalizer localizers[3][3];
  PlanarSeparator separators[3][3];
  LocalizedSeparatorLink link_combined[3][3];
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      cell[j][i] = unit_cell;
      cell[j][i].shift(static_cast<double>(i - 1), static_cast<double>(j - 1),
                       0.0);
      localizers[j][i] = cell[j][i].getLocalizer();
    }
  }
  Normal plane_normals[3][3];
  plane_normals[0][0] = Normal(-1.0, -1.0, 0.0);
  plane_normals[0][2] = Normal(1.0, -1.0, 0.0);
  plane_normals[2][0] = Normal(-1.0, 1.0, 0.0);
  plane_normals[2][2] = Normal(1.0, 1.0, 0.0);
  plane_normals[0][1] = Normal(0.0, -1.0, 0.0);
  plane_normals[2][1] = Normal(0.0, 1.0, 0.0);
  plane_normals[1][0] = Normal(-1.0, 0.0, 0.0);
  plane_normals[1][2] = Normal(1.0, 0.0, 0.0);
  plane_normals[1][1] = Normal(0.0, 0.0, 0.0);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      plane_normals[j][i].normalize();
      plane_normals[j][i] = -plane_normals[j][i];
      localizers[j][i].addPlane(
          Plane(plane_normals[j][i],
                plane_normals[j][i] * (cell[j][i].calculateCentroid())));
      separators[j][i] =
          PlanarSeparator::fromOnePlane(Plane(Normal(-1.0, 0.0, 0.0), 0.0));
      link_combined[j][i] =
          LocalizedSeparatorLink(&localizers[j][i], &separators[j][i]);
    }
  }
  (localizers[1][1])[6].distance() = (-1000000.00);

  SeparatedMoments<VolumeMoments> correct_volume_moments[3][3];
  double ss = 7.0 / 6.0;
  correct_volume_moments[0][0] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.50, Pt(-ss, -ss, 0.0)));
  correct_volume_moments[0][1] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.25, Pt(0.25, -1.25, 0.0)),
      VolumeMoments(0.25, Pt(-0.25, -1.25, 0.0)));
  correct_volume_moments[0][2] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.50, Pt(ss, -ss, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));
  correct_volume_moments[1][0] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.50, Pt(-1.25, 0.0, 0.0)));
  correct_volume_moments[1][1] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));
  correct_volume_moments[1][2] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.50, Pt(1.25, 0.0, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));
  correct_volume_moments[2][0] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.50, Pt(-ss, ss, 0.0)));
  correct_volume_moments[2][1] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.25, Pt(0.25, 1.25, 0.0)),
      VolumeMoments(0.25, Pt(-0.25, 1.25, 0.0)));
  correct_volume_moments[2][2] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.50, Pt(ss, ss, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));

  // Check individual cutting by these LocalizedSeparators
  SeparatedMoments<VolumeMoments> individual_volume_moments[9];
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      individual_volume_moments[j * 3 + i] =
          getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
              super_cell, link_combined[j][i]);
    }
  }
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      UnsignedIndex_t ind = static_cast<UnsignedIndex_t>(j * 3 + i);
      EXPECT_NEAR(individual_volume_moments[ind][0].volume(),
                  correct_volume_moments[j][i][0].volume(), 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].volume(),
                  correct_volume_moments[j][i][1].volume(), 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][0].centroid()[0],
                  correct_volume_moments[j][i][0].centroid()[0], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][0].centroid()[1],
                  correct_volume_moments[j][i][0].centroid()[1], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][0].centroid()[2],
                  correct_volume_moments[j][i][0].centroid()[2], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].centroid()[0],
                  correct_volume_moments[j][i][1].centroid()[0], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].centroid()[1],
                  correct_volume_moments[j][i][1].centroid()[1], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].centroid()[2],
                  correct_volume_moments[j][i][1].centroid()[2], 1.0e-14);
    }
  }

  LocalizedSeparatorLink* ptr;
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      link_combined[j][i].setId(static_cast<UnsignedIndex_t>(j * 3 + i));
      ptr = i - 1 >= 0 ? &link_combined[j][i - 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(0, ptr);
      ptr = i + 1 < 3 ? &link_combined[j][i + 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(1, ptr);
      ptr = j - 1 >= 0 ? &link_combined[j - 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(2, ptr);
      ptr = j + 1 < 3 ? &link_combined[j + 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(3, ptr);
    }
  }

  // Should be agnostic to starting location
  auto volume_moments = getNormalizedVolumeMoments<
      AccumulatedVolumeMoments<SeparatedMoments<VolumeMoments>>>(
      super_cell, link_combined[0][0]);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      auto ind = static_cast<UnsignedIndex_t>(j * 3 + i);
    }
  }
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      auto ind = static_cast<UnsignedIndex_t>(j * 3 + i);
      EXPECT_NEAR(volume_moments[ind][0].volume(),
                  correct_volume_moments[j][i][0].volume(), 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].volume(),
                  correct_volume_moments[j][i][1].volume(), 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[0],
                  correct_volume_moments[j][i][0].centroid()[0], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[1],
                  correct_volume_moments[j][i][0].centroid()[1], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[2],
                  correct_volume_moments[j][i][0].centroid()[2], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[0],
                  correct_volume_moments[j][i][1].centroid()[0], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[1],
                  correct_volume_moments[j][i][1].centroid()[1], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[2],
                  correct_volume_moments[j][i][1].centroid()[2], 1.0e-14);
    }
  }

  // Should be agnostic to starting location
  volume_moments = getNormalizedVolumeMoments<
      AccumulatedVolumeMoments<SeparatedMoments<VolumeMoments>>>(
      super_cell, link_combined[2][2]);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      auto ind = static_cast<UnsignedIndex_t>(j * 3 + i);
      EXPECT_NEAR(volume_moments[ind][0].volume(),
                  correct_volume_moments[j][i][0].volume(), 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].volume(),
                  correct_volume_moments[j][i][1].volume(), 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[0],
                  correct_volume_moments[j][i][0].centroid()[0], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[1],
                  correct_volume_moments[j][i][0].centroid()[1], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][0].centroid()[2],
                  correct_volume_moments[j][i][0].centroid()[2], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[0],
                  correct_volume_moments[j][i][1].centroid()[0], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[1],
                  correct_volume_moments[j][i][1].centroid()[1], 1.0e-14);
      EXPECT_NEAR(volume_moments[ind][1].centroid()[2],
                  correct_volume_moments[j][i][1].centroid()[2], 1.0e-14);
    }
  }
}

TEST(GenericCutting, MostlyAlignedVolume) {
  Dodecahedron flux_volume(
      {Pt(-2.0000000000000000, 0.0000000000000000, -2.0000000000000000),
       Pt(-2.0000000000000000, -0.50000000000000000, -2.0000000000000000),
       Pt(-1.5000000000000000, -0.50000000000000000, -2.0000000000000000),
       Pt(-1.5000000000000000, 0.0000000000000000, -2.0000000000000000),
       Pt(-2.0000000000000000, 0.0000000000000000, -1.9875000000000000),
       Pt(-2.0000000000000000, -0.50000000000000000, -1.9875000000000000),
       Pt(-1.5000000000000000, -0.50000000000000000, -1.9875000000000000),
       Pt(-1.5000000000000000, 0.0000000000000000, -1.9875000000000000)});
  auto flux_moments = flux_volume.calculateMoments();
  flux_moments.normalizeByVolume();

  RectangularCuboid enclosing_cell = RectangularCuboid::fromBoundingPts(
      Pt(-2.0, -0.5, -2.0), Pt(-1.5, 0.0, -1.5));

  PlanarLocalizer localizer_for_cell = enclosing_cell.getLocalizer();
  PlanarSeparator below_everything_separator =
      PlanarSeparator::fromOnePlane(Plane(Normal(1.0, 0.0, 0.0), 1000000.0));
  LocalizedSeparatorLink localized_separator(&localizer_for_cell,
                                             &below_everything_separator);

  VolumeMoments separated_moments =
      getNormalizedVolumeMoments<VolumeMoments, RecursiveSimplexCutting>(
          flux_volume, localizer_for_cell);
}

TEST(GenericCutting, PlanarSeparatorPathTriangle) {
  static constexpr UnsignedIndex_t nlayers = 1;
  // static constexpr Pt bottom_left_vertex = Pt(3.22, 4.69, 0.0);
  static constexpr Pt bottom_left_vertex = Pt(-1.0, 0.0, 0.0);
  static constexpr double width = 2.0;
  static constexpr double height = 1.0;
  std::array<std::array<Pt, 3>, 2> layer_points;
  const auto top_point =
      Pt(bottom_left_vertex[0] + 0.5 * width, bottom_left_vertex[1] + height,
         bottom_left_vertex[2]);

  auto slanted_normal =
      Normal::fromPtNormalized(top_point - bottom_left_vertex);
  auto flipped_normal = slanted_normal;
  flipped_normal[0] = -flipped_normal[0];
  const auto upwards_pointing_normal = Normal(0.0, 1.0, 0.0);

  // Initialize base points
  layer_points[0][0] = bottom_left_vertex;
  layer_points[0][2] = bottom_left_vertex + Pt(width, 0.0, 0.0);
  layer_points[0][1] = 0.5 * (layer_points[0][0] + layer_points[0][2]);

  auto base_triangle = Tri({layer_points[0][0], layer_points[0][2], top_point});
  base_triangle.setPlaneOfExistence(
      Plane(Normal(0.0, 0.0, 1.0), bottom_left_vertex[2]));

  std::vector<Tri> triangles;
  std::vector<PlanarSeparator> separators;
  triangles.reserve(3 * nlayers + 1);
  separators.reserve(3 * nlayers + 1);
  for (UnsignedIndex_t layer = 0; layer < nlayers; ++layer) {
    layer_points[1][0] = 0.5 * (layer_points[0][0] + top_point);
    layer_points[1][2] = 0.5 * (layer_points[0][2] + top_point);
    layer_points[1][1] = 0.5 * (layer_points[1][0] + layer_points[1][2]);

    // Add triangles
    triangles.push_back(
        Tri({layer_points[0][0], layer_points[0][1], layer_points[1][0]}));
    triangles.push_back(
        Tri({layer_points[0][1], layer_points[0][2], layer_points[1][2]}));
    triangles.push_back(
        Tri({layer_points[0][1], layer_points[1][2], layer_points[1][0]}));

    // Add Planes that will create this volume
    separators.push_back(PlanarSeparator::fromOnePlane(
        Plane(slanted_normal, slanted_normal * layer_points[0][1])));
    separators.push_back(PlanarSeparator::fromOnePlane(
        Plane(flipped_normal, flipped_normal * layer_points[0][1])));
    separators.push_back(PlanarSeparator::fromOnePlane(
        Plane(upwards_pointing_normal,
              upwards_pointing_normal * layer_points[1][1])));

    layer_points[0] = layer_points[1];
  }

  // Add left over top triangle and encompassing plane
  triangles.push_back(Tri({layer_points[0][0], layer_points[0][2], top_point}));
  separators.push_back(
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 0.0), 1.0)));

  std::vector<PlanarSeparatorPath> separator_path(separators.size());
  for (UnsignedIndex_t n = 0; n < separators.size(); ++n) {
    separator_path[n] = PlanarSeparatorPath(&separators[n]);
    separator_path[n].setId(n);
  }
  for (UnsignedIndex_t n = 0; n < separators.size() - 1; ++n) {
    separator_path[n].setEdgeConnectivity(&separator_path[n + 1]);
  }

  auto moments =
      getNormalizedVolumeMoments<TaggedAccumulatedVolumeMoments<VolumeMoments>>(
          base_triangle, separator_path[0]);
  for (UnsignedIndex_t n = 0; n < triangles.size(); ++n) {
    // Should be direct correspondence between n, the uniqueId, and the
    // triangle;
    triangles[n].setPlaneOfExistence(
        Plane(Normal(0.0, 0.0, 1.0), bottom_left_vertex[2]));
    auto triangle_moments = triangles[n].calculateMoments();
    triangle_moments.normalizeByVolume();
    EXPECT_NEAR(moments[n].volume(), triangle_moments.volume(), 1.0e-14);
    EXPECT_NEAR(moments[n].centroid()[0], triangle_moments.centroid()[0],
                1.0e-14);
    EXPECT_NEAR(moments[n].centroid()[1], triangle_moments.centroid()[1],
                1.0e-14);
    EXPECT_NEAR(moments[n].centroid()[2], triangle_moments.centroid()[2],
                1.0e-14);
  }

  PlanarSeparatorPathGroup path_group;
  for (UnsignedIndex_t n = 0; n < separators.size(); ++n) {
    path_group.addPlanarSeparatorPath(separator_path[n]);
  }
  std::vector<UnsignedIndex_t> priority(3 * nlayers + 1);
  std::iota(priority.begin(), priority.end(), 0);
  path_group.setPriorityOrder(priority);

  moments =
      getNormalizedVolumeMoments<TaggedAccumulatedVolumeMoments<VolumeMoments>>(
          base_triangle, path_group.getFirstReconstruction());
  for (UnsignedIndex_t n = 0; n < triangles.size(); ++n) {
    // Should be direct correspondence between n, the uniqueId, and the
    // triangle;
    triangles[n].setPlaneOfExistence(
        Plane(Normal(0.0, 0.0, 1.0), bottom_left_vertex[2]));
    auto triangle_moments = triangles[n].calculateMoments();
    triangle_moments.normalizeByVolume();
    EXPECT_NEAR(moments[n].volume(), triangle_moments.volume(), 1.0e-14);
    EXPECT_NEAR(moments[n].centroid()[0], triangle_moments.centroid()[0],
                1.0e-14);
    EXPECT_NEAR(moments[n].centroid()[1], triangle_moments.centroid()[1],
                1.0e-14);
    EXPECT_NEAR(moments[n].centroid()[2], triangle_moments.centroid()[2],
                1.0e-14);
  }
}

TEST(GenericCutting, getVolumeMomentsPlanarSeparatorPathGroupTagged) {
  RectangularCuboid super_cell = RectangularCuboid::fromBoundingPts(
      Pt(-1.5, -1.5, -0.5), Pt(1.5, 1.5, 0.5));
  RectangularCuboid cell[3][3];
  PlanarLocalizer localizers[3][3];
  PlanarSeparator separators[3][3][2];
  PlanarSeparatorPathGroup separator_group[3][3];
  LocalizedSeparatorGroupLink link_combined[3][3];
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      cell[j][i] = unit_cell;
      cell[j][i].shift(static_cast<double>(i - 1), static_cast<double>(j - 1),
                       0.0);
      localizers[j][i] = cell[j][i].getLocalizer();
    }
  }
  Normal plane_normals[3][3];
  plane_normals[0][0] = Normal(-1.0, -1.0, 0.0);
  plane_normals[0][2] = Normal(1.0, -1.0, 0.0);
  plane_normals[2][0] = Normal(-1.0, 1.0, 0.0);
  plane_normals[2][2] = Normal(1.0, 1.0, 0.0);
  plane_normals[0][1] = Normal(0.0, -1.0, 0.0);
  plane_normals[2][1] = Normal(0.0, 1.0, 0.0);
  plane_normals[1][0] = Normal(-1.0, 0.0, 0.0);
  plane_normals[1][2] = Normal(1.0, 0.0, 0.0);
  plane_normals[1][1] = Normal(0.0, 0.0, 0.0);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      plane_normals[j][i].normalize();
      plane_normals[j][i] = -plane_normals[j][i];
      localizers[j][i].addPlane(
          Plane(plane_normals[j][i],
                plane_normals[j][i] * (cell[j][i].calculateCentroid())));
      separators[j][i][0] =
          PlanarSeparator::fromOnePlane(Plane(Normal(-1.0, 0.0, 0.0), 0.0));
      separators[j][i][1] =
          PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 0.0), 1.0));

      separator_group[j][i].addPlanarSeparatorPath(
          PlanarSeparatorPath(&separators[j][i][0]), 0);
      separator_group[j][i].addPlanarSeparatorPath(
          PlanarSeparatorPath(&separators[j][i][1]), 1);
      separator_group[j][i].setPriorityOrder({0, 1});

      link_combined[j][i] = LocalizedSeparatorGroupLink(&localizers[j][i],
                                                        &separator_group[j][i]);
    }
  }
  (localizers[1][1])[6].distance() = (-1000000.00);

  SeparatedMoments<VolumeMoments> correct_volume_moments[3][3];
  double ss = 7.0 / 6.0;
  correct_volume_moments[0][0] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.50, Pt(-ss, -ss, 0.0)));
  correct_volume_moments[0][1] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.25, Pt(0.25, -1.25, 0.0)),
      VolumeMoments(0.25, Pt(-0.25, -1.25, 0.0)));
  correct_volume_moments[0][2] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.50, Pt(ss, -ss, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));
  correct_volume_moments[1][0] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.50, Pt(-1.25, 0.0, 0.0)));
  correct_volume_moments[1][1] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));
  correct_volume_moments[1][2] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.50, Pt(1.25, 0.0, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));
  correct_volume_moments[2][0] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)),
                                      VolumeMoments(0.50, Pt(-ss, ss, 0.0)));
  correct_volume_moments[2][1] = SeparatedMoments<VolumeMoments>(
      VolumeMoments(0.25, Pt(0.25, 1.25, 0.0)),
      VolumeMoments(0.25, Pt(-0.25, 1.25, 0.0)));
  correct_volume_moments[2][2] =
      SeparatedMoments<VolumeMoments>(VolumeMoments(0.50, Pt(ss, ss, 0.0)),
                                      VolumeMoments(0.00, Pt(0.0, 0.0, 0.0)));

  // Check individual cutting by these LocalizedSeparatorGroup
  AccumulateWrapper<TaggedAccumulatedVolumeMoments<VolumeMoments>>
      individual_volume_moments[9];
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      individual_volume_moments[j * 3 + i] = getNormalizedVolumeMoments<
          AccumulateWrapper<TaggedAccumulatedVolumeMoments<VolumeMoments>>>(
          super_cell, link_combined[j][i]);
    }
  }
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      auto ind = static_cast<UnsignedIndex_t>(j * 3 + i);
      EXPECT_NEAR(individual_volume_moments[ind][0].volume(),
                  correct_volume_moments[j][i][0].volume(), 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].volume(),
                  correct_volume_moments[j][i][1].volume(), 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][0].centroid()[0],
                  correct_volume_moments[j][i][0].centroid()[0], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][0].centroid()[1],
                  correct_volume_moments[j][i][0].centroid()[1], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][0].centroid()[2],
                  correct_volume_moments[j][i][0].centroid()[2], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].centroid()[0],
                  correct_volume_moments[j][i][1].centroid()[0], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].centroid()[1],
                  correct_volume_moments[j][i][1].centroid()[1], 1.0e-14);
      EXPECT_NEAR(individual_volume_moments[ind][1].centroid()[2],
                  correct_volume_moments[j][i][1].centroid()[2], 1.0e-14);
    }
  }

  LocalizedSeparatorGroupLink* ptr;
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      link_combined[j][i].setId(
          static_cast<UnsignedIndex_t>(8 * j * 3 + i));  // Non-contigous tag
      ptr = i - 1 >= 0 ? &link_combined[j][i - 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(0, ptr);
      ptr = i + 1 < 3 ? &link_combined[j][i + 1] : nullptr;
      link_combined[j][i].setEdgeConnectivity(1, ptr);
      ptr = j - 1 >= 0 ? &link_combined[j - 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(2, ptr);
      ptr = j + 1 < 3 ? &link_combined[j + 1][i] : nullptr;
      link_combined[j][i].setEdgeConnectivity(3, ptr);
    }
  }

  // Should be agnostic to starting location
  auto volume_moments =
      getNormalizedVolumeMoments<TaggedAccumulatedVolumeMoments<
          TaggedAccumulatedVolumeMoments<VolumeMoments>>>(super_cell,
                                                          link_combined[0][0]);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      // Lookup via known tag value
      auto tag = static_cast<UnsignedIndex_t>(8 * j * 3 + i);
      EXPECT_NEAR(volume_moments[tag][0].volume(),
                  correct_volume_moments[j][i][0].volume(), 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][1].volume(),
                  correct_volume_moments[j][i][1].volume(), 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][0].centroid()[0],
                  correct_volume_moments[j][i][0].centroid()[0], 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][0].centroid()[1],
                  correct_volume_moments[j][i][0].centroid()[1], 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][0].centroid()[2],
                  correct_volume_moments[j][i][0].centroid()[2], 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][1].centroid()[0],
                  correct_volume_moments[j][i][1].centroid()[0], 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][1].centroid()[1],
                  correct_volume_moments[j][i][1].centroid()[1], 1.0e-14);
      EXPECT_NEAR(volume_moments[tag][1].centroid()[2],
                  correct_volume_moments[j][i][1].centroid()[2], 1.0e-14);
    }
  }

  // Should be agnostic to starting location
  volume_moments = getNormalizedVolumeMoments<TaggedAccumulatedVolumeMoments<
      TaggedAccumulatedVolumeMoments<VolumeMoments>>>(super_cell,
                                                      link_combined[2][2]);
  for (UnsignedIndex_t elem = 0; elem < volume_moments.size(); ++elem) {
    auto tag = volume_moments.getTagForIndex(elem);
    int j = static_cast<int>(tag / 24);
    int i = static_cast<int>(tag) - j * 24;
    // Get tag value for entry and compare to correct moments
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[0].volume(),
                correct_volume_moments[j][i][0].volume(), 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[1].volume(),
                correct_volume_moments[j][i][1].volume(), 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[0].centroid()[0],
                correct_volume_moments[j][i][0].centroid()[0], 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[0].centroid()[1],
                correct_volume_moments[j][i][0].centroid()[1], 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[0].centroid()[2],
                correct_volume_moments[j][i][0].centroid()[2], 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[1].centroid()[0],
                correct_volume_moments[j][i][1].centroid()[0], 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[1].centroid()[1],
                correct_volume_moments[j][i][1].centroid()[1], 1.0e-14);
    EXPECT_NEAR(volume_moments.getMomentsForIndex(elem)[1].centroid()[2],
                correct_volume_moments[j][i][1].centroid()[2], 1.0e-14);
  }
}

TEST(GenericCutting, SmallShearedVolume) {
  Dodecahedron dodecahedron({Pt(15.0, 7.0, 4.0), Pt(15.0, 8.0, 4.0),
                             Pt(15.0, 8.0, 5.0), Pt(15.0, 7.0, 5.0),
                             Pt(15.0 + 0.1 * DBL_EPSILON, 7.0, 4.5),
                             Pt(15.0 + 0.1 * DBL_EPSILON, 8.0, 4.5),
                             Pt(15.0 + 0.1 * DBL_EPSILON, 8.0, 4.75),
                             Pt(15.0 + 0.1 * DBL_EPSILON, 7.0, 4.75)});
  auto cell = RectangularCuboid::fromBoundingPts(Pt(15.0, 7.0, 4.0),
                                                 Pt(16.5, 8.0, 5.0));
  auto localizer = cell.getLocalizer();
  auto left_cell = RectangularCuboid::fromBoundingPts(Pt(14.0, 7.0, 4.0),
                                                      Pt(15.0, 8.0, 5.0));
  auto left_localizer = left_cell.getLocalizer();
  LocalizerLink localizer_link(&localizer);
  LocalizerLink left_localizer_link(&left_localizer);
  localizer_link.setEdgeConnectivity(0, &left_localizer_link);
  localizer_link.setId(0);
  left_localizer_link.setEdgeConnectivity(1, &localizer_link);
  left_localizer_link.setId(1);

  Volume vol = getNormalizedVolumeMoments<Volume, HalfEdgeCutting>(
      dodecahedron, localizer_link);
  Volume from_left_vol = getNormalizedVolumeMoments<Volume, HalfEdgeCutting>(
      dodecahedron, left_localizer_link);
}

// TEST(GenericCutting, MichaelsProblem) {
//  CappedDodecahedron flux_volume(
//      {Pt(0.10937500000000000, 0.32812500000000000, 0.15625000000000000),
//       Pt(9.3750000000000000E-002, 0.32812500000000000, 0.15625000000000000),
//       Pt(9.3750000000000000E-002, 0.32812500000000000, 0.17187500000000000),
//       Pt(0.10937500000000000, 0.32812500000000000, 0.17187500000000000),
//       Pt(0.11011970283132823, 0.32817864557724924, 0.15640127741437027),
//       Pt(9.4322656794659343E-002, 0.32810287829862494, 0.15595822134122375),
//       Pt(9.3608558426119287E-002, 0.32807167591881881, 0.17109713574351781),
//       Pt(0.10873015159778411, 0.32822855987369370, 0.17215871525186258),
//       Pt(0.10244680569966691, 0.32809663029964958, 0.16381904516875029)});
//
//  RectangularCuboid enclosing_cell = RectangularCuboid::fromBoundingPts(
//      Pt(0.09375, 0.3125, 0.140625), Pt(0.109375, 0.328125, 0.15625));
//  PlanarLocalizer localizer_for_cell = enclosing_cell.getLocalizer();
//  PlanarSeparator below_everything_separator =
//      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 0.0), -1.0));
//  LocalizedSeparatorLink localized_separator(&localizer_for_cell,
//                                             &below_everything_separator);
//  setMinimumVolumeToTrack(1.0e-15 * enclosing_cell.calculateVolume());
//
//  auto correct_moments = flux_volume.calculateMoments();
//  correct_moments.normalizeByVolume();
//  std::cout << "Correct moments: " << correct_moments << std::endl;
//
//  //  auto cut_moments =
//  //  getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>,
//  //                                                HalfEdgeCutting>(
//  //      flux_volume, localized_separator);
//  auto cut_moments = getVolumeMoments<VolumeMoments, HalfEdgeCutting>(
//      flux_volume, localizer_for_cell);
//
//  // auto &complete_polytope = setHalfEdgeStructure(flux_volume);
//  //  auto half_edge_polytope =
//  //      generateSegmentedVersion<CappedDodecahedron>(&complete_polytope);
//  //  assert(half_edge_polytope.checkValidHalfEdgeStructure());
//  //
//  //  auto cut_moments = getVolumeMoments<VolumeMoments, HalfEdgeCutting>(
//  //      &half_edge_polytope, &complete_polytope, localizer_for_cell);
//  //
//  //  std::cout << half_edge_polytope << std::endl;
//
//  std::cout << "UnNormalized Cut moments: " << cut_moments << std::endl;
//  cut_moments.normalizeByVolume();
//  std::cout << "Cut moments: " << cut_moments << std::endl;
//}

TEST(GenericCutting, CircularHanging) {
  Dodecahedron flux_volume =
      Dodecahedron({Pt(-2.38612715975352, -2.00005071922666, -2.30840410376319),
                    Pt(-2.37459844293369, -3.0, -2.31872130759989),
                    Pt(-1.99745667999418, -3.0, -3.0),
                    Pt(-2.01554412994774, -1.98191255004644, -3.0),
                    Pt(-2.39612715975352, -2.00005071922666, -2.30840410376319),
                    Pt(-2.38459844293369, -3.0, -2.31872130759989),
                    Pt(-2.00745667999418, -3.0, -3.0),
                    Pt(-2.02554412994774, -1.98191255004644, -3.0)});

  PlanarLocalizer localizers[5];
  PlanarSeparator separators[5];
  LocalizedSeparatorLink localizer_link[5];

  localizers[0].setNumberOfPlanes(5);
  localizers[0][0] =
      Plane(Normal(0.878089656813471, 0.010440881766227, 0.478382213909697),
            -3.22566800969335);
  localizers[0][1] =
      Plane(Normal(-0.742287107880983, -0.00772563686508022, 0.670037434781606),
            0.23992592989337);
  localizers[0][2] = Plane(Normal(0.0, 0.0, -1.0), 3.0);
  localizers[0][3] =
      Plane(Normal(0.0157851068272069, 0.99971765082297, 0.0177608849851179),
            -2.00282003421817);
  localizers[0][4] = Plane(Normal(-0.0, -1.0, 0.0), 3.0);
  localizer_link[0] = LocalizedSeparatorLink(&localizers[0], &separators[0]);
  localizer_link[0].setId(968);
  localizer_link[0].setEdgeConnectivity(0, &localizer_link[1]);

  localizers[1].setNumberOfPlanes(5);
  localizers[1][0] =
      Plane(Normal(-0.321804003032829, -0.0178264802260627, 0.94663847388261),
            -1.38170500153929);
  localizers[1][1] =
      Plane(Normal(-0.878089656813471, -0.010440881766227, -0.478382213909697),
            3.22566800969335);
  localizers[1][2] =
      Plane(Normal(0.930283406450889, 0.0220313997169217, -0.36617946571134),
            -0.820153169519917);
  localizers[1][3] =
      Plane(Normal(-0.0157842188348016, 0.999717656387519, 0.0177613609507049),
            -2.00282326293126);
  localizers[1][4] = Plane(Normal(-0.0, -1.0, 0.0), 3.0);
  localizer_link[1] = LocalizedSeparatorLink(&localizers[1], &separators[1]);
  localizer_link[1].setId(971);
  localizer_link[1].setEdgeConnectivity(0, &localizer_link[2]);

  localizers[2].setNumberOfPlanes(5);
  localizers[2][0] =
      Plane(Normal(0.737357139929398, -0.000319989284423859, 0.675503031674911),
            -2.59426543441096);
  localizers[2][1] =
      Plane(Normal(-0.986768615447812, -0.0150505436821806, 0.161434756796913),
            2.01200058934869);
  localizers[2][2] =
      Plane(Normal(0.321804003032829, 0.0178264802260627, -0.94663847388261),
            1.38170500153929);
  localizers[2][3] =
      Plane(Normal(-0.015784365891251, 0.999717646517383, 0.0177617858098666),
            -2.00282387304173);
  localizers[2][4] = Plane(Normal(-0.0, -1.0, 0.0), 3.0);
  localizer_link[2] = LocalizedSeparatorLink(&localizers[2], &separators[2]);
  localizer_link[2].setId(972);
  localizer_link[2].setEdgeConnectivity(1, &localizer_link[3]);

  localizers[3].setNumberOfPlanes(5);
  localizers[3][0] =
      Plane(Normal(-0.602371152393203, -0.0205693529865265, 0.797951061458157),
            0.267102672821628);
  localizers[3][1] =
      Plane(Normal(-0.463296977615825, 0.00885707713853614, -0.886158824769354),
            3.13338356550321);
  localizers[3][2] =
      Plane(Normal(0.986768615447812, 0.0150505436821806, -0.161434756796913),
            -2.01200058934869);
  localizers[3][3] =
      Plane(Normal(-0.0157849200594005, 0.999717636179301, 0.0177618752047237),
            -2.00282273640882);
  localizers[3][4] = Plane(Normal(-0.0, -1.0, 0.0), 3.0);
  localizer_link[3] = LocalizedSeparatorLink(&localizers[3], &separators[3]);
  localizer_link[3].setId(982);
  localizer_link[3].setEdgeConnectivity(1, &localizer_link[4]);

  localizers[4].setNumberOfPlanes(5);
  localizers[4][0] =
      Plane(Normal(0.463296977615825, -0.00885707713853614, 0.886158824769354),
            -3.13338356550321);
  localizers[4][1] = Plane(Normal(-1.0, 0.0, 0.0), 3.0);
  localizers[4][2] =
      Plane(Normal(0.742287107880983, 0.00772563686508022, -0.670037434781606),
            -0.23992592989337);
  localizers[4][3] =
      Plane(Normal(-0.0157853378506084, 0.999717643532565, 0.0177610900182606),
            -2.00281994168512);
  localizers[4][4] = Plane(Normal(-0.0, -1.0, 0.0), 3.0);
  localizer_link[4] = LocalizedSeparatorLink(&localizers[4], &separators[4]);
  localizer_link[4].setId(969);
  localizer_link[4].setEdgeConnectivity(2, &localizer_link[0]);

  Volume test_vol = getNormalizedVolumeMoments<Volume, HalfEdgeCutting>(
      flux_volume, localizer_link[0]);
}

TEST(GenericCutting, MethodDifference) {
  auto cap_octa = CappedOctahedron_LLL(
      {Pt(0.0, 0.40625, -0.015625), Pt(0.03125, 0.375, 0.015625),
       Pt(0.0, 0.375, 0.015625),
       Pt(0.0252686335297825, 0.405993193219365, -0.015625),
       Pt(0.0549982680028043, 0.373222827692387, 0.015625),
       Pt(0.0237482680028043, 0.375186323100881, 0.015625),
       Pt(0.0345813716481551, 0.383362792071091, 0.00378633922404048)});

  PlanarLocalizer localizer;
  localizer.addPlane(Plane(Normal(0.0, 0.707106781186547, 0.707106781186547),
                           0.276213586400995));
  localizer.addPlane(Plane(Normal(-0.707106781186547, -0.0, 0.707106781186547),
                           -0.0110485434560398));
  localizer.addPlane(
      Plane(Normal(0.577350269189626, -0.577350269189626, -0.577350269189626),
            -0.207485252990022));
  localizer.addPlane(Plane(Normal(0.0, 0.0, -1.0), 0.015625));

  PlanarSeparator separator;
  separator.addPlane(
      Plane(Normal(-0.207018186962485, -0.976966875955128, -0.0517609268970713),
            -0.392491275508266));

  LocalizedSeparator localized_separator(&localizer, &separator);

  auto simplex_moments =
      getVolumeMoments<SeparatedMoments<Volume>, SimplexCutting>(
          cap_octa, localized_separator);
  auto recursive_moments =
      getVolumeMoments<SeparatedMoments<Volume>, RecursiveSimplexCutting>(
          cap_octa, localized_separator);

  EXPECT_NEAR(simplex_moments[0], recursive_moments[0], 1.0e-14);
  EXPECT_NEAR(simplex_moments[1], recursive_moments[1], 1.0e-14);
}
}  // namespace
