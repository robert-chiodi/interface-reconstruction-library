// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/planar_reconstruction/planar_separator.h"
#include "gtest/gtest.h"

#include "src/generic_cutting/generic_cutting.h"
#include "src/geometry/general/plane.h"
#include "src/interface_reconstruction_methods/reconstruction_cleaning.h"

namespace {

using namespace IRL;

TEST(Separators, SeparatorExtra) {
  // Construction and flipping of cut phase
  PlanarSeparator reconstruction;
  EXPECT_EQ(reconstruction.isFlipped(), false);
  reconstruction.flipCutting();
  EXPECT_EQ(reconstruction.isFlipped(), true);

  // Construction of two plane reconstruction
  Plane plane_1(Normal::normalized(1.0, 1.0, 1.0), 0.25);
  Plane plane_2(Normal::normalized(-1.0, -1.0, -1.0), 0.25);
  PlanarSeparator recon_two_plane =
      PlanarSeparator::fromTwoPlanes(plane_1, plane_2, 1.0);
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

  // Removing all planes and setting distance large to cause single phase cell
  plane_1 = Plane(Normal::normalized(1.0, 1.0, 1.0), 0.25);
  plane_2 = Plane(Normal::normalized(-1.0, -1.0, -1.0), 0.25);
  recon_two_plane = PlanarSeparator::fromTwoPlanes(plane_1, plane_2, 1.0);
  recon_two_plane.flipCutting();
  setToPurePhaseReconstruction(1.0, &recon_two_plane);
  EXPECT_EQ(recon_two_plane.isFlipped(), false);
  EXPECT_EQ(recon_two_plane.getNumberOfPlanes(), 1);
  Plane zeroed_plane = recon_two_plane[0];
  EXPECT_EQ(zeroed_plane.distance() > 0.0, true);
  setToPurePhaseReconstruction(0.0, &recon_two_plane);
  zeroed_plane = recon_two_plane[0];
  EXPECT_EQ(zeroed_plane.distance() < 0.0, true);

  // Ability to clean remove and add planes
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

  PlanarSeparator flip_setting =
      PlanarSeparator::fromTwoPlanes(Plane(Normal(1.0, 0.0, 0.0), 0.25),
                                     Plane(Normal(0.0, 1.0, 0.0), 0.25), 1.0);
  flip_setting.setFlip(-1.0);
  EXPECT_EQ(flip_setting.isFlipped(), true);

  // Test serialization
  ByteBuffer buffer;
  std::array<PlanarSeparator, 2> planar_reconstructions_to_pack;
  planar_reconstructions_to_pack[0].addPlane(
      Plane(Normal::normalized(4.8, -2.2, 6.9), -13.08));
  planar_reconstructions_to_pack[0].addPlane(
      Plane(Normal::normalized(-4.8, 2.2, -6.9), -13.08));
  planar_reconstructions_to_pack[0].flipCutting();
  PlanarLocalizer cube_localizer = unit_cell.getLocalizer();
  for (UnsignedIndex_t n = 0; n < 6; ++n) {
    planar_reconstructions_to_pack[1].addPlane(cube_localizer[n]);
  }
  serializeAndPack(planar_reconstructions_to_pack[0], &buffer);
  serializeAndPack(planar_reconstructions_to_pack[1], &buffer);
  buffer.resetBufferPointer();
  std::array<PlanarSeparator, 2> recv_planar_reconstruction;
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
