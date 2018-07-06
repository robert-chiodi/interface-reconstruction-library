// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/polygons/tri.h"

#include <random>

#include "gtest/gtest.h"

#include "src/generic_cutting/cut_polygon.h"
#include "src/geometry/general/normal.h"
#include "src/geometry/general/rotations.h"
#include "src/geometry/polygons/divided_polygon.h"
#include "src/helpers/helper.h"

namespace {

using namespace IRL;

TEST(Tri, Tri) {
  Pt pt_list[3] = {Pt(-1.0, 0.0, 0.0), Pt(1.0, 0.0, 0.0), Pt(0.0, 1.0, 0.0)};
  Tri tri_3pt({pt_list[0], pt_list[1], pt_list[2]});
  EXPECT_DOUBLE_EQ(tri_3pt[0].x(), -1.0);
  EXPECT_DOUBLE_EQ(tri_3pt[0].y(), 0.0);
  EXPECT_DOUBLE_EQ(tri_3pt[0].z(), 0.0);
  EXPECT_DOUBLE_EQ(tri_3pt[1].x(), 1.0);
  EXPECT_DOUBLE_EQ(tri_3pt[1].y(), 0.0);
  EXPECT_DOUBLE_EQ(tri_3pt[1].z(), 0.0);
  EXPECT_DOUBLE_EQ(tri_3pt[2].x(), 0.0);
  EXPECT_DOUBLE_EQ(tri_3pt[2].y(), 1.0);
  EXPECT_DOUBLE_EQ(tri_3pt[2].z(), 0.0);

  Tri tri_list = Tri::fromRawPtPointer(3, pt_list);
  for (UnsignedIndex_t v = 0; v < 3; ++v) {
    EXPECT_DOUBLE_EQ(tri_list[v].x(), tri_3pt[v].x());
    EXPECT_DOUBLE_EQ(tri_list[v].y(), tri_3pt[v].y());
    EXPECT_DOUBLE_EQ(tri_list[v].z(), tri_3pt[v].z());
  }

  Normal tri_normal = tri_3pt.calculateNormal();
  EXPECT_DOUBLE_EQ(tri_normal[0], 0.0);
  EXPECT_DOUBLE_EQ(tri_normal[1], 0.0);
  EXPECT_DOUBLE_EQ(tri_normal[2], 1.0);

  tri_3pt.calculateAndSetPlaneOfExistence();
  Pt tri_centroid = tri_3pt.calculateCentroid();
  EXPECT_DOUBLE_EQ(tri_centroid[0], 0.0);
  EXPECT_DOUBLE_EQ(tri_centroid[1], 1.0 / 3.0);
  EXPECT_DOUBLE_EQ(tri_centroid[2], 0.0);

  // Test serialization
  ByteBuffer buffer;
  std::array<Tri, 5> tri_to_pack;
  DividedPolygon polygon;
  polygon.addVertex(Pt(0.0, 0.0, 0.0));
  polygon.addVertex(Pt(1.0, 0.0, 0.0));
  polygon.addVertex(Pt(1.0, 1.0, 0.0));
  polygon.addVertex(Pt(0.0, 1.0, 0.0));
  polygon.addVertex(Pt(-1.0, 1.0, 0.0));
  polygon.calculateAndSetPlaneOfExistence();
  polygon.resetCentroid();

  Pt centroid = polygon.calculateCentroid();

  tri_to_pack[0] = Tri({polygon[0], centroid, polygon[4]});
  tri_to_pack[1] = Tri({polygon[1], centroid, polygon[0]});
  tri_to_pack[2] = Tri({polygon[2], centroid, polygon[1]});
  tri_to_pack[3] = Tri({polygon[3], centroid, polygon[2]});
  tri_to_pack[4] = Tri({polygon[4], centroid, polygon[3]});
  serializeAndPack(tri_to_pack[0], &buffer);
  serializeAndPack(tri_to_pack[1], &buffer);
  serializeAndPack(tri_to_pack[2], &buffer);
  serializeAndPack(tri_to_pack[3], &buffer);
  serializeAndPack(tri_to_pack[4], &buffer);
  buffer.resetBufferPointer();
  std::array<Tri, 5> recv_tris;
  unpackAndStore(&recv_tris[0], &buffer);
  unpackAndStore(&recv_tris[1], &buffer);
  unpackAndStore(&recv_tris[2], &buffer);
  unpackAndStore(&recv_tris[3], &buffer);
  unpackAndStore(&recv_tris[4], &buffer);
  for (UnsignedIndex_t n = 0; n < 5; ++n) {
    for (UnsignedIndex_t v = 0; v < 3; ++v) {
      EXPECT_DOUBLE_EQ(recv_tris[n][v].x(), tri_to_pack[n][v].x());
      EXPECT_DOUBLE_EQ(recv_tris[n][v].y(), tri_to_pack[n][v].y());
      EXPECT_DOUBLE_EQ(recv_tris[n][v].z(), tri_to_pack[n][v].z());
    }
  }
}

}  // namespace
