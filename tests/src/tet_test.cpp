// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/polyhedrons/tet.h"

#include "gtest/gtest.h"

#include <algorithm>

#include "src/geometry/general/pt.h"
#include "src/geometry/general/rotations.h"
#include "src/geometry/polygons/tri.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/helpers/geometric_cutting_helpers.h"

namespace {

using namespace IRL;

TEST(Polyhedra, Tet) {
  // Construction of tet and access
  Tet set_tet({Pt(-1.0, -1.0, -1.0), Pt(-1.0, 0.0, 0.0), Pt(0.0, 0.0, -1.0),
               Pt(0.0, 0.0, 0.0)});
  EXPECT_DOUBLE_EQ(set_tet[0].x(), -1.0);
  EXPECT_DOUBLE_EQ(set_tet[0].y(), -1.0);
  EXPECT_DOUBLE_EQ(set_tet[0].z(), -1.0);
  EXPECT_DOUBLE_EQ(set_tet[1].x(), -1.0);
  EXPECT_DOUBLE_EQ(set_tet[1].y(), 0.0);
  EXPECT_DOUBLE_EQ(set_tet[1].z(), 0.0);
  EXPECT_DOUBLE_EQ(set_tet[2].x(), 0.0);
  EXPECT_DOUBLE_EQ(set_tet[2].y(), 0.0);
  EXPECT_DOUBLE_EQ(set_tet[2].z(), -1.0);
  EXPECT_DOUBLE_EQ(set_tet[3].x(), 0.0);
  EXPECT_DOUBLE_EQ(set_tet[3].y(), 0.0);
  EXPECT_DOUBLE_EQ(set_tet[3].z(), 0.0);

  // Calculation of geometric quantities
  double vol = set_tet.calculateVolume();
  Pt centroid = set_tet.calculateCentroid();
  EXPECT_DOUBLE_EQ(vol, 1.0 / 6.0);
  EXPECT_DOUBLE_EQ(centroid[0], -0.5);
  EXPECT_DOUBLE_EQ(centroid[1], -0.25);
  EXPECT_DOUBLE_EQ(centroid[2], -0.5);

  VolumeMoments tet_moments = set_tet.calculateMoments();
  EXPECT_DOUBLE_EQ(tet_moments.volume(), vol);
  EXPECT_DOUBLE_EQ(tet_moments.centroid()[0], vol * centroid[0]);
  EXPECT_DOUBLE_EQ(tet_moments.centroid()[1], vol * centroid[1]);
  EXPECT_DOUBLE_EQ(tet_moments.centroid()[2], vol * centroid[2]);

  double flat_pt_list[12] = {-1.0, -1.0, -1.0, -1.0, 0.0, 0.0,
                             0.0,  0.0,  -1.0, 0.0,  0.0, 0.0};
  Tet list_set_tet = Tet::fromRawDoublePointer(4, flat_pt_list);
  VolumeMoments pt_list_moments = list_set_tet.calculateMoments();
  EXPECT_DOUBLE_EQ(pt_list_moments.volume(), vol);
  EXPECT_DOUBLE_EQ(pt_list_moments.centroid()[0], vol * centroid[0]);
  EXPECT_DOUBLE_EQ(pt_list_moments.centroid()[1], vol * centroid[1]);
  EXPECT_DOUBLE_EQ(pt_list_moments.centroid()[2], vol * centroid[2]);

  Tet negative_tet({Pt(-0.5, -0.5, 0.5), Pt(-0.5, 0.5, 0.5),
                    Pt(-0.5, 0.5, -0.5), Pt(-0.8, 0.5, 0.5)});
  EXPECT_DOUBLE_EQ(negative_tet.calculateSign(), -1.0);

  Tet positive_tet = negative_tet;
  positive_tet[3] = Pt(0.0, 0.5, 0.5);
  EXPECT_DOUBLE_EQ(positive_tet.calculateSign(), 1.0);

  // Test serialization
  RectangularCuboid two_point_set_rect_cubic =
      RectangularCuboid::fromBoundingPts(Pt(-0.5, -0.5, -0.5),
                                         Pt(0.5, 0.5, 0.5));
  ByteBuffer buffer;
  std::array<Tet, 5> tets_to_pack;
  tets_to_pack[0] = unit_cell.getSimplexFromDecomposition(0);
  tets_to_pack[1] = unit_cell.getSimplexFromDecomposition(1);
  tets_to_pack[2] = unit_cell.getSimplexFromDecomposition(2);
  tets_to_pack[3] = unit_cell.getSimplexFromDecomposition(3);
  tets_to_pack[4] = unit_cell.getSimplexFromDecomposition(4);
  serializeAndPack(tets_to_pack[0], &buffer);
  serializeAndPack(tets_to_pack[1], &buffer);
  serializeAndPack(tets_to_pack[2], &buffer);
  serializeAndPack(tets_to_pack[3], &buffer);
  serializeAndPack(tets_to_pack[4], &buffer);
  buffer.resetBufferPointer();
  std::array<Tet, 5> recv_tets;
  unpackAndStore(&recv_tets[0], &buffer);
  unpackAndStore(&recv_tets[1], &buffer);
  unpackAndStore(&recv_tets[2], &buffer);
  unpackAndStore(&recv_tets[3], &buffer);
  unpackAndStore(&recv_tets[4], &buffer);
  for (UnsignedIndex_t n = 0; n < 5; ++n) {
    for (UnsignedIndex_t v = 0; v < 4; ++v) {
      EXPECT_DOUBLE_EQ(recv_tets[n][v].x(), tets_to_pack[n][v].x());
      EXPECT_DOUBLE_EQ(recv_tets[n][v].y(), tets_to_pack[n][v].y());
      EXPECT_DOUBLE_EQ(recv_tets[n][v].z(), tets_to_pack[n][v].z());
    }
  }
}
}  // namespace
