// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/half_edge_structures/half_edge_polytope.h"

#include "gtest/gtest.h"

#include "src/generic_cutting/half_edge_cutting/half_edge_cutting.h"
#include "src/geometry/polyhedrons/concave_box.h"
#include "src/geometry/polyhedrons/dodecahedron.h"
#include "src/geometry/polyhedrons/hexahedron.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/geometry/polyhedrons/tet.h"
#include "src/helpers/geometric_cutting_helpers.h"
#include "src/moments/volume.h"

namespace {

using namespace IRL;

TEST(HalfEdgePolytope, Debug) {
  Tet tet_to_write({Pt(1.0, 0.0, -0.5), Pt(1.0, 0.0, 0.5), Pt(1.0, 1.0, 0.0),
                    Pt(0.0, 0.5, 0.0)});
  auto tet_half_edge = tet_to_write.generateHalfEdgeVersion();
  auto tet_segmented_half_edge = tet_half_edge.generateSegmentedPolyhedron();
  assert(tet_segmented_half_edge.checkValidHalfEdgeStructure());

  Volume tmp_volume;
  PlanarLocalizer cutting_reconstruction =
      PlanarLocalizer::fromOnePlane(Plane(Normal(1.0, 0.0, 0.0), 0.5));
  decltype(tet_segmented_half_edge) clipped_tet;
  splitHalfEdgePolytope(&tet_segmented_half_edge, &clipped_tet, &tet_half_edge,
                        cutting_reconstruction[0]);
  assert(tet_segmented_half_edge.checkValidHalfEdgeStructure());

  Dodecahedron rect_cuboid;
  for (UnsignedIndex_t n = 0; n < 8; ++n) {
    rect_cuboid[n] = unit_cell[n];
  }
  rect_cuboid[6][0] = 1.0;
  rect_cuboid[7][0] = 1.0;
  rect_cuboid[4][0] = 0.0;
  rect_cuboid[5][0] = 0.0;
  auto rc_half_edge = rect_cuboid.generateHalfEdgeVersion();
  auto rc_segmented_half_edge = rc_half_edge.generateSegmentedPolyhedron();
  assert(rc_segmented_half_edge.checkValidHalfEdgeStructure());
  cutting_reconstruction[0] = Plane(Normal::normalized(1.0, 1.0, 0.0), 0.0);
  decltype(rc_segmented_half_edge) clipped_polyhedron;
  splitHalfEdgePolytope(&rc_segmented_half_edge, &clipped_polyhedron,
                        &rc_half_edge, cutting_reconstruction[0]);
  Dodecahedron pointy;
  for (UnsignedIndex_t n = 0; n < 8; ++n) {
    pointy[n] = unit_cell[n];
  }
  pointy[2][1] = 1.0;
  pointy[7][1] = 1.0;
  auto pointy_half_edge = pointy.generateHalfEdgeVersion();
  auto pointy_segmented_half_edge =
      pointy_half_edge.generateSegmentedPolyhedron();
  assert(pointy_segmented_half_edge.checkValidHalfEdgeStructure());
  cutting_reconstruction[0] = Plane(Normal::normalized(0.0, 1.0, 0.0), 0.6);
  decltype(pointy_segmented_half_edge) clipped_pointy;

  splitHalfEdgePolytope(&pointy_segmented_half_edge, &clipped_pointy,
                        &pointy_half_edge, cutting_reconstruction[0]);
  ConcaveBox conbox;
  conbox[0] = Pt(1.0, 0.0, 0.0);
  conbox[1] = Pt(1.0, 1.0, 0.0);
  conbox[2] = Pt(1.0, 1.0, 0.25);
  conbox[3] = Pt(0.5, 1.0, 0.3);
  conbox[4] = Pt(1.0, 1.0, 0.35);
  conbox[5] = Pt(1.0, 1.0, 0.6);
  conbox[6] = Pt(1.0, 0.0, 0.6);
  conbox[7] = Pt(1.0, 0.0, 0.35);
  conbox[8] = Pt(0.5, 0.0, 0.3);
  conbox[9] = Pt(1.0, 0.0, 0.25);
  conbox[10] = Pt(0.0, 0.0, 0.0);
  conbox[11] = Pt(0.0, 1.0, 0.0);
  conbox[12] = Pt(0.0, 1.0, 0.6);
  conbox[13] = Pt(0.0, 0.0, 0.6);
  auto conbox_volume_moments = conbox.calculateMoments();
  conbox_volume_moments.normalizeByVolume();

  EXPECT_NEAR(conbox_volume_moments.volume(), 0.575, 1.0e-14);
  EXPECT_NEAR(conbox_volume_moments.centroid()[0], 67.0 / 138.0, 1.0e-14);
  EXPECT_NEAR(conbox_volume_moments.centroid()[1], 0.5, 1.0e-14);
  EXPECT_NEAR(conbox_volume_moments.centroid()[2], 0.3, 1.0e-14);
  auto conv_half_edge = conbox.generateHalfEdgeVersion();
  auto conv_segmented = conv_half_edge.generateSegmentedPolyhedron();
  assert(conv_segmented.checkValidHalfEdgeStructure());
  cutting_reconstruction[0] = Plane(Normal::normalized(1.0, 0.0, 0.0), 0.75);
  decltype(conv_segmented) clipped_convbox;
  splitHalfEdgePolytope(&conv_segmented, &clipped_convbox, &conv_half_edge,
                        cutting_reconstruction[0]);
  assert(conv_segmented.checkValidHalfEdgeStructure());
  assert(clipped_convbox.checkValidHalfEdgeStructure());
  auto half_edge_volume_moments = conv_segmented.calculateMoments();
  half_edge_volume_moments.normalizeByVolume();
  EXPECT_NEAR(half_edge_volume_moments.volume(), 0.44375, 1.0e-14);
  EXPECT_NEAR(half_edge_volume_moments.centroid()[0], 79.0 / 213.0, 1.0e-14);
  EXPECT_NEAR(half_edge_volume_moments.centroid()[1], 0.5, 1.0e-14);
  EXPECT_NEAR(half_edge_volume_moments.centroid()[2], 0.3, 1.0e-14);
  half_edge_volume_moments = clipped_convbox.calculateMoments();
  half_edge_volume_moments.normalizeByVolume();
  EXPECT_NEAR(half_edge_volume_moments.volume(), 0.13125, 1.0e-14);
  EXPECT_NEAR(half_edge_volume_moments.centroid()[0], 55.0 / 63.0, 1.0e-14);
  EXPECT_NEAR(half_edge_volume_moments.centroid()[1], 0.5, 1.0e-14);
  EXPECT_NEAR(half_edge_volume_moments.centroid()[2], 0.3, 1.0e-14);
}

TEST(HalfEdgePolygon, Debug) {
  Polygon pentagon;
  pentagon.addVertex(Pt(0.0, 0.0, 0.0));
  pentagon.addVertex(Pt(1.0, 0.0, 0.0));
  pentagon.addVertex(Pt(1.0, 1.0, 0.0));
  pentagon.addVertex(Pt(0.5, 1.5, 0.0));
  pentagon.addVertex(Pt(0.0, 1.0, 0.0));
  auto pentagon_half_edge = pentagon.generateHalfEdgeVersion();
  auto segmented_pentagon_half_edge =
      pentagon_half_edge.generateSegmentedPolygon();
  PlanarLocalizer cutting_reconstruction = PlanarLocalizer::fromOnePlane(
      Plane(Normal::normalized(1.0, 0.0, 0.0), 0.5));
  decltype(segmented_pentagon_half_edge) clipped_penta;
  splitHalfEdgePolytope(&segmented_pentagon_half_edge, &clipped_penta,
                        &pentagon_half_edge, cutting_reconstruction[0]);

  Polygon crossed_quad;
  crossed_quad.addVertex(Pt(0.0, 0.0, 0.0));
  crossed_quad.addVertex(Pt(1.0, 0.0, 0.0));
  crossed_quad.addVertex(Pt(0.0, 1.0, 0.0));
  crossed_quad.addVertex(Pt(1.0, 1.0, 0.0));
  auto cq_half_edge = crossed_quad.generateHalfEdgeVersion();
  auto segmented_cq_half_edge = cq_half_edge.generateSegmentedPolygon();
  cutting_reconstruction = PlanarLocalizer::fromOnePlane(
      Plane(Normal::normalized(1.0, 0.0, 0.0), 0.5));
  decltype(segmented_cq_half_edge) clipped_cq;
  splitHalfEdgePolytope(&segmented_cq_half_edge, &clipped_cq, &cq_half_edge,
                        cutting_reconstruction[0]);
}
}  // namespace
