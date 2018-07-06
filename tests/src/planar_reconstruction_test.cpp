// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/planar_reconstruction/planar_reconstruction.h"
#include "gtest/gtest.h"

#include "src/geometry/general/normal.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"

namespace {

using namespace IRL;
TEST(PlanarReconstruction, PlanarReconstruction) {
  // Construction of one plane reconstruction
  Plane plane_1(Normal::normalized(1.0, 1.0, 1.0), 0.25);
  Plane plane_2(Normal::normalized(-1.0, -1.0, -1.0), 0.25);
  PlanarReconstruction<1> recon_one_plane =
      PlanarReconstruction<1>::fromOnePlane(plane_1);
  EXPECT_EQ(recon_one_plane.getNumberOfPlanes(), 1);
  Plane plane_1_out = recon_one_plane[0];
  EXPECT_DOUBLE_EQ(plane_1_out.normal()[0], plane_1.normal()[0]);
  EXPECT_DOUBLE_EQ(plane_1_out.normal()[1], plane_1.normal()[1]);
  EXPECT_DOUBLE_EQ(plane_1_out.normal()[2], plane_1.normal()[2]);
  EXPECT_DOUBLE_EQ(plane_1_out.distance(), plane_1.distance());

  // Construction of two plane reconstruction
  PlanarReconstruction<2> recon_two_plane =
      PlanarReconstruction<2>::fromTwoPlanes(plane_1, plane_2);
  EXPECT_EQ(recon_two_plane.getNumberOfPlanes(), 2);
  Plane plane_2_out = recon_two_plane[1];
  EXPECT_DOUBLE_EQ(plane_2_out.normal()[0], plane_2.normal()[0]);
  EXPECT_DOUBLE_EQ(plane_2_out.normal()[1], plane_2.normal()[1]);
  EXPECT_DOUBLE_EQ(plane_2_out.normal()[2], plane_2.normal()[2]);
  EXPECT_DOUBLE_EQ(plane_2_out.distance(), plane_2.distance());
  EXPECT_DOUBLE_EQ(recon_two_plane[1].normal()[0], plane_2.normal()[0]);
  EXPECT_DOUBLE_EQ(recon_two_plane[1].normal()[1], plane_2.normal()[1]);
  EXPECT_DOUBLE_EQ(recon_two_plane[1].normal()[2], plane_2.normal()[2]);
  EXPECT_DOUBLE_EQ(recon_two_plane[1].distance(), plane_2.distance());

  // Resetting a plane
  recon_two_plane[1] = plane_1;
  plane_1_out = recon_two_plane[0];
  plane_2_out = recon_two_plane[1];
  EXPECT_DOUBLE_EQ(plane_1_out.normal()[0], plane_2_out.normal()[0]);
  EXPECT_DOUBLE_EQ(plane_1_out.normal()[1], plane_2_out.normal()[1]);
  EXPECT_DOUBLE_EQ(plane_1_out.normal()[2], plane_2_out.normal()[2]);
  EXPECT_DOUBLE_EQ(plane_1_out.distance(), plane_2_out.distance());

  // Setting just distance
  SmallVector<double, 2> dist;
  dist.resize(2);
  dist[0] = 3.2;
  dist[1] = 6.8;
  recon_two_plane.setDistances(dist);
  EXPECT_DOUBLE_EQ(recon_two_plane[0].distance(), dist[0]);
  EXPECT_DOUBLE_EQ(recon_two_plane[1].distance(), dist[1]);
  recon_two_plane[0].distance() = 8.6;
  EXPECT_DOUBLE_EQ(recon_two_plane[0].distance(), 8.6);

  // Resetting of number of interfaces
  recon_two_plane.zeroNumberOfPlanes();
  EXPECT_EQ(recon_two_plane.getNumberOfPlanes(), 0);

  // Test serialization
  ByteBuffer buffer;
  std::array<PlanarReconstruction<2>, 2> planar_reconstructions_to_pack;
  planar_reconstructions_to_pack[0].addPlane(
      Plane(Normal::normalized(4.8, -2.2, 6.9), 13.08));
  planar_reconstructions_to_pack[0].addPlane(
      Plane(Normal::normalized(-4.8, 2.2, -6.9), 13.08));
  planar_reconstructions_to_pack[1].setNumberOfPlanes(6);
  for (UnsignedIndex_t n = 0; n < 6; ++n) {
    planar_reconstructions_to_pack[1][n] = unit_cell.getLocalizer()[n];
  }

  serializeAndPack(planar_reconstructions_to_pack[0], &buffer);
  serializeAndPack(planar_reconstructions_to_pack[1], &buffer);
  buffer.resetBufferPointer();
  std::array<PlanarReconstruction<2>, 2> recv_planar_reconstruction;
  unpackAndStore(&recv_planar_reconstruction[0], &buffer);
  unpackAndStore(&recv_planar_reconstruction[1], &buffer);
  for (UnsignedIndex_t n = 0; n < 2; ++n) {
    EXPECT_EQ(recv_planar_reconstruction[n].getNumberOfPlanes(),
              planar_reconstructions_to_pack[n].getNumberOfPlanes());
    for (UnsignedIndex_t p = 0;
         p < planar_reconstructions_to_pack[n].getNumberOfPlanes(); ++p) {
      EXPECT_DOUBLE_EQ(recv_planar_reconstruction[n][p].normal()[0],
                       planar_reconstructions_to_pack[n][p].normal()[0]);
      EXPECT_DOUBLE_EQ(recv_planar_reconstruction[n][p].normal()[1],
                       planar_reconstructions_to_pack[n][p].normal()[1]);
      EXPECT_DOUBLE_EQ(recv_planar_reconstruction[n][p].normal()[2],
                       planar_reconstructions_to_pack[n][p].normal()[2]);
      EXPECT_DOUBLE_EQ(recv_planar_reconstruction[n][p].distance(),
                       planar_reconstructions_to_pack[n][p].distance());
    }
  }
}

}  // namespace
