// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/interface_reconstruction_methods/reconstruction_interface.h"
#include "src/planar_reconstruction/planar_separator.h"

#include <cmath>
#include <random>

#include <fstream>
#include <iomanip>

#include "gtest/gtest.h"
#include "src/generic_cutting/generic_cutting.h"

#include "src/generic_cutting/cut_polygon.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/general/rotations.h"
#include "src/helpers/mymath.h"
#include "src/interface_reconstruction_methods/plane_distance.h"

namespace {

using namespace IRL;

TEST(ReconstructionInterface, R2P_2D1P) {
  RectangularCuboid middle_cell = RectangularCuboid::fromBoundingPts(
      Pt(90.0, 90.0, 90.0), Pt(100.0, 100.0, 100.0));
  Normal plane_normal = Normal::normalized(-1.0, 1.0, 0.0);
  Plane correct_plane(
      plane_normal,
      plane_normal * Pt(middle_cell.calculateCentroid() - Pt(2.0, 1.0, 0.0)));
  PlanarSeparator correct_reconstruction =
      PlanarSeparator::fromOnePlane(correct_plane);
  R2PNeighborhood<RectangularCuboid> neighborhood_geometry;
  neighborhood_geometry.resize(9);
  RectangularCuboid cells[9];
  SeparatedMoments<VolumeMoments> svm[9];

  for (int j = -1; j < 2; ++j) {
    for (int i = -1; i < 2; ++i) {
      UnsignedIndex_t ind = static_cast<UnsignedIndex_t>((j + 1) * 3 + i + 1);
      cells[ind] = middle_cell;
      cells[ind].shift(10.0 * static_cast<double>(i),
                       10.0 * static_cast<double>(j), 0.0);
      svm[ind] = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          cells[ind], correct_reconstruction);
      neighborhood_geometry.setMember(ind, &cells[ind], &svm[ind]);
    }
  }
  neighborhood_geometry.setSurfaceArea(
      getReconstructionSurfaceArea(middle_cell, correct_reconstruction));
  neighborhood_geometry.setCenterOfStencil(4);

  // Perturb initial reconstruction then see if it is recovered.
  Plane initial_plane(Normal(-1.0, 0.0, 0.0),
                      Normal(-1.0, 0.0, 0.0) * middle_cell.calculateCentroid());
  PlanarSeparator initial_reconstruction =
      PlanarSeparator::fromOnePlane(initial_plane);
  setDistanceToMatchVolumeFractionPartialFill(
      neighborhood_geometry.getCenterCell(),
      (neighborhood_geometry.getCenterCellStoredMoments())[0].volume() /
          neighborhood_geometry.getCenterCell().calculateVolume(),
      &initial_reconstruction);
  PlanarSeparator final_reconstruction =
      reconstructionWithR2P2D(neighborhood_geometry, initial_reconstruction);
  double rotation_angle;
  Normal rotation_axis;
  rotateNormalOntoNormal(final_reconstruction[0].normal(),
                         correct_plane.normal(), &rotation_angle,
                         &rotation_axis);

  EXPECT_NEAR(getVolumeFraction(middle_cell, final_reconstruction),
              (neighborhood_geometry.getCenterCellStoredMoments())[0].volume() /
                  neighborhood_geometry.getCenterCell().calculateVolume(),
              1.0e-14);

  std::cout
      << "\n\n"
         "******** Correctness is not automatically checked here ********\n"
      << "Target normal: " << correct_reconstruction[0].normal()[0] << " "
      << correct_reconstruction[0].normal()[1] << " "
      << correct_reconstruction[0].normal()[2] << '\n'
      << "Target distance: " << correct_reconstruction[0].distance() << '\n'
      << "Found normal:  " << final_reconstruction[0].normal()[0] << " "
      << final_reconstruction[0].normal()[1] << " "
      << final_reconstruction[0].normal()[2] << '\n'
      << "Found distance: " << final_reconstruction[0].distance() << '\n'
      << "The angle between the two of these is: " << rad2Deg(rotation_angle)
      << " degrees \n"
      << "******** Check this ^^^^^^^^ for correctness ********\n"
      << "To print out path of optimization (in rotated frame)\n"
      << "change reconstructionWithR2P2D to reconstructionWithR2P2DDebug\n\n";

  // Same problem on unit cell below this

  middle_cell = unit_cell;
  plane_normal = Normal(-1.0, 1.0, 0.0);
  plane_normal.normalize();
  correct_plane = Plane(
      plane_normal,
      plane_normal * Pt(middle_cell.calculateCentroid() - Pt(0.2, 0.1, 0.0)));
  correct_reconstruction = PlanarSeparator::fromOnePlane(correct_plane);

  for (int j = -1; j < 2; ++j) {
    for (int i = -1; i < 2; ++i) {
      UnsignedIndex_t ind = static_cast<UnsignedIndex_t>((j + 1) * 3 + i + 1);
      cells[ind] = middle_cell;
      cells[ind].shift(static_cast<double>(i), static_cast<double>(j), 0.0);
      svm[ind] = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          cells[ind], correct_reconstruction);
      neighborhood_geometry.setMember(ind, &cells[ind], &svm[ind]);
    }
  }
  neighborhood_geometry.setSurfaceArea(
      getReconstructionSurfaceArea(middle_cell, correct_reconstruction));
  neighborhood_geometry.setCenterOfStencil(4);

  // Perturb initial reconstruction then see if it is recovered.
  initial_plane =
      Plane(Normal(-1.0, 0.0, 0.0),
            Normal(-1.0, 0.0, 0.0) * middle_cell.calculateCentroid());
  initial_reconstruction = PlanarSeparator::fromOnePlane(initial_plane);
  setDistanceToMatchVolumeFractionPartialFill(
      neighborhood_geometry.getCenterCell(),
      (neighborhood_geometry.getCenterCellStoredMoments())[0].volume() /
          neighborhood_geometry.getCenterCell().calculateVolume(),
      &initial_reconstruction);
  final_reconstruction =
      reconstructionWithR2P2D(neighborhood_geometry, initial_reconstruction);
  rotateNormalOntoNormal(final_reconstruction[0].normal(),
                         correct_plane.normal(), &rotation_angle,
                         &rotation_axis);

  EXPECT_NEAR(getVolumeFraction(middle_cell, final_reconstruction),
              (neighborhood_geometry.getCenterCellStoredMoments())[0].volume() /
                  neighborhood_geometry.getCenterCell().calculateVolume(),
              1.0e-14);

  std::cout
      << "\n\n"
         "******** Correctness is not automatically checked here ********\n"
      << "Target normal: " << correct_reconstruction[0].normal()[0] << " "
      << correct_reconstruction[0].normal()[1] << " "
      << correct_reconstruction[0].normal()[2] << '\n'
      << "Target distance: " << correct_reconstruction[0].distance() << '\n'
      << "Found normal:  " << final_reconstruction[0].normal()[0] << " "
      << final_reconstruction[0].normal()[1] << " "
      << final_reconstruction[0].normal()[2] << '\n'
      << "Found distance: " << final_reconstruction[0].distance() << '\n'
      << "The angle between the two of these is: " << rad2Deg(rotation_angle)
      << " degrees \n"
      << "******** Check this ^^^^^^^^ for correctness ********\n"
      << "To print out path of optimization (in rotated frame)\n"
      << "change reconstructionWithR2P2D to reconstructionWithR2P2DDebug\n\n";
}

TEST(ReconstructionInterface, R2P_3D1P) {
  Plane correct_plane(Normal::normalized(-1.0, 1.0, -1.0), 0.25);
  PlanarSeparator correct_reconstruction =
      PlanarSeparator::fromOnePlane(correct_plane);
  R2PNeighborhood<RectangularCuboid> neighborhood_geometry;
  neighborhood_geometry.resize(27);
  RectangularCuboid cells[27];
  SeparatedMoments<VolumeMoments> svm[27];
  for (int k = -1; k < 2; ++k) {
    for (int j = -1; j < 2; ++j) {
      for (int i = -1; i < 2; ++i) {
        UnsignedIndex_t ind =
            static_cast<UnsignedIndex_t>((k + 1) * 9 + (j + 1) * 3 + i + 1);
        cells[ind] = unit_cell;
        cells[ind].shift(static_cast<double>(i), static_cast<double>(j),
                         static_cast<double>(k));
        svm[ind] = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
            cells[ind], correct_reconstruction);
        neighborhood_geometry.setMember(ind, &cells[ind], &svm[ind]);
      }
    }
  }
  neighborhood_geometry.setSurfaceArea(
      getReconstructionSurfaceArea(unit_cell, correct_reconstruction));
  neighborhood_geometry.setCenterOfStencil(13);

  // Perturb initial reconstruction then see if it is recovered.
  Plane initial_plane(Normal(-1.0, 0.0, 0.0), 0.0);
  PlanarSeparator initial_reconstruction =
      PlanarSeparator::fromOnePlane(initial_plane);
  setDistanceToMatchVolumeFractionPartialFill(
      unit_cell,
      (neighborhood_geometry.getCenterCellStoredMoments())[0].volume(),
      &initial_reconstruction);
  PlanarSeparator final_reconstruction;
  final_reconstruction =
      reconstructionWithR2P3D(neighborhood_geometry, initial_reconstruction);
  double rotation_angle;
  Normal rotation_axis;
  rotateNormalOntoNormal(final_reconstruction[0].normal(),
                         correct_plane.normal(), &rotation_angle,
                         &rotation_axis);

  EXPECT_NEAR(getVolumeFraction(unit_cell, final_reconstruction),
              (neighborhood_geometry.getCenterCellStoredMoments())[0].volume(),
              1.0e-14);

  std::cout
      << "\n\n"
         "******** Correctness is not automatically checked here ********\n"
      << "Target normal: " << correct_reconstruction[0].normal()[0] << " "
      << correct_reconstruction[0].normal()[1] << " "
      << correct_reconstruction[0].normal()[2] << '\n'
      << "Found normal:  " << final_reconstruction[0].normal()[0] << " "
      << final_reconstruction[0].normal()[1] << " "
      << final_reconstruction[0].normal()[2] << '\n'
      << "The angle between the two of these is: " << rad2Deg(rotation_angle)
      << " degrees \n"
      << "******** Check this ^^^^^^^^ for correctness ********\n"
      << "To print out path of optimization (in rotated frame)\n"
      << "change reconstructionWithR2P2D to reconstructionWithR2P3DDebug\n\n";
}

TEST(ReconstructionInterface, R2P_2D2P) {
  RectangularCuboid middle_cell = RectangularCuboid::fromBoundingPts(
      Pt(-100.0, -100.0, -100.0), Pt(-80.0, -80.0, -80.0));
  Normal normal_0 = Normal::normalized(-1.0, 1.0, 0.0);
  Normal normal_1 = Normal::normalized(1.0, 1.0, 0.0);
  Plane correct_plane_0(normal_0, normal_0 * middle_cell.calculateCentroid());
  Plane correct_plane_1(normal_1, normal_1 * middle_cell.calculateCentroid());
  PlanarSeparator correct_reconstruction =
      PlanarSeparator::fromTwoPlanes(correct_plane_0, correct_plane_1, -1.0);
  R2PNeighborhood<RectangularCuboid> neighborhood_geometry;
  neighborhood_geometry.resize(9);
  RectangularCuboid cells[9];
  SeparatedMoments<VolumeMoments> svm[9];
  for (int j = -1; j < 2; ++j) {
    for (int i = -1; i < 2; ++i) {
      UnsignedIndex_t ind = static_cast<UnsignedIndex_t>((j + 1) * 3 + i + 1);
      cells[ind] = middle_cell;
      cells[ind].shift(static_cast<double>(i) * 20.0,
                       static_cast<double>(j) * 20.0, 0.0);
      svm[ind] = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          cells[ind], correct_reconstruction);
      neighborhood_geometry.setMember(ind, &cells[ind], &svm[ind]);
    }
  }
  neighborhood_geometry.setSurfaceArea(
      getReconstructionSurfaceArea(middle_cell, correct_reconstruction));
  neighborhood_geometry.setCenterOfStencil(4);

  // Perturb initial reconstruction then see if it is recovered.
  Normal rotation_axis = Normal::normalized(0.0, 0.0, 1.0);
  double rotation_angle = deg2Rad(45.0);
  Normal initial_normal_0 =
      UnitQuaternion(rotation_angle, rotation_axis) * correct_plane_0.normal();
  Normal initial_normal_1 =
      UnitQuaternion(rotation_angle, rotation_axis) * correct_plane_1.normal();
  Plane initial_plane_0(initial_normal_0,
                        initial_normal_0 * middle_cell.calculateCentroid());
  Plane initial_plane_1(initial_normal_1,
                        initial_normal_1 * middle_cell.calculateCentroid());
  PlanarSeparator initial_reconstruction =
      PlanarSeparator::fromTwoPlanes(initial_plane_0, initial_plane_1, -1.0);
  setDistanceToMatchVolumeFractionPartialFill(
      middle_cell,
      (neighborhood_geometry.getCenterCellStoredMoments())[0].volume() /
          middle_cell.calculateVolume(),
      &initial_reconstruction);
  PlanarSeparator final_reconstruction;
  final_reconstruction =
      reconstructionWithR2P2D(neighborhood_geometry, initial_reconstruction);
  rotateNormalOntoNormal(final_reconstruction[0].normal(),
                         correct_plane_0.normal(), &rotation_angle,
                         &rotation_axis);

  EXPECT_EQ(final_reconstruction.isFlipped(),
            correct_reconstruction.isFlipped());

  EXPECT_EQ(final_reconstruction.getNumberOfPlanes(), 2);

  EXPECT_NEAR(getVolumeFraction(middle_cell, final_reconstruction),
              (neighborhood_geometry.getCenterCellStoredMoments())[0].volume() /
                  middle_cell.calculateVolume(),
              global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);

  std::cout
      << "\n\n"
         "******** Correctness is not automatically checked here ********\n"
      << "Target normal[0]: " << correct_reconstruction[0].normal()[0] << " "
      << correct_reconstruction[0].normal()[1] << " "
      << correct_reconstruction[0].normal()[2] << '\n'
      << "Target normal[1]: " << correct_reconstruction[1].normal()[0] << " "
      << correct_reconstruction[1].normal()[1] << " "
      << correct_reconstruction[1].normal()[2] << '\n'
      << "Target distances: " << correct_reconstruction[0].distance() << " "
      << correct_reconstruction[1].distance() << '\n'
      << "Initial normal[0]: " << initial_reconstruction[0].normal()[0] << " "
      << initial_reconstruction[0].normal()[1] << " "
      << initial_reconstruction[0].normal()[2] << '\n'
      << "Initial normal[1]: " << initial_reconstruction[1].normal()[0] << " "
      << initial_reconstruction[1].normal()[1] << " "
      << initial_reconstruction[1].normal()[2] << '\n'
      << "Initial distances: " << initial_reconstruction[0].distance() << " "
      << initial_reconstruction[1].distance() << '\n'
      << "Found normal[0]:  " << final_reconstruction[0].normal()[0] << " "
      << final_reconstruction[0].normal()[1] << " "
      << final_reconstruction[0].normal()[2] << '\n'
      << "Found normal[1]:  " << final_reconstruction[1].normal()[0] << " "
      << final_reconstruction[1].normal()[1] << " "
      << final_reconstruction[1].normal()[2] << '\n'
      << "Found distances: " << final_reconstruction[0].distance() << " "
      << final_reconstruction[1].distance() << '\n'
      << "The angle between the two shared normals is: "
      << rad2Deg(rotation_angle) << " degrees \n"
      << "******** Check this ^^^^^^^^ for correctness ********\n"
      << "To print out path of optimization (in rotated frame)\n"
      << "change reconstructionWithR2P2D to reconstructionWithR2P2DDebug\n\n";

  // Unit cell variant of above test below.
  middle_cell = unit_cell;
  normal_0 = Normal::normalized(-1.0, 1.0, 0.0);
  normal_1 = Normal::normalized(1.0, 1.0, 0.0);
  correct_plane_0 = Plane(normal_0, normal_0 * middle_cell.calculateCentroid());
  correct_plane_1 = Plane(normal_1, normal_1 * middle_cell.calculateCentroid());
  correct_reconstruction =
      PlanarSeparator::fromTwoPlanes(correct_plane_0, correct_plane_1, -1.0);

  for (int j = -1; j < 2; ++j) {
    for (int i = -1; i < 2; ++i) {
      UnsignedIndex_t ind = static_cast<UnsignedIndex_t>((j + 1) * 3 + i + 1);
      cells[ind] = middle_cell;
      cells[ind].shift(static_cast<double>(i), static_cast<double>(j), 0.0);
      svm[ind] = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
          cells[ind], correct_reconstruction);
      neighborhood_geometry.setMember(ind, &cells[ind], &svm[ind]);
    }
  }
  neighborhood_geometry.setSurfaceArea(
      getReconstructionSurfaceArea(middle_cell, correct_reconstruction));
  neighborhood_geometry.setCenterOfStencil(4);

  // Perturb initial reconstruction then see if it is recovered.
  rotation_axis = Normal(0.0, 0.0, 1.0);
  rotation_angle = deg2Rad(45.0);
  initial_normal_0 =
      UnitQuaternion(rotation_angle, rotation_axis) * correct_plane_0.normal();
  initial_normal_1 =
      UnitQuaternion(rotation_angle, rotation_axis) * correct_plane_1.normal();
  initial_plane_0 = Plane(initial_normal_0,
                          initial_normal_0 * middle_cell.calculateCentroid());
  initial_plane_1 = Plane(initial_normal_1,
                          initial_normal_1 * middle_cell.calculateCentroid());
  initial_reconstruction =
      PlanarSeparator::fromTwoPlanes(initial_plane_0, initial_plane_1, -1.0);
  setDistanceToMatchVolumeFractionPartialFill(
      middle_cell,
      (neighborhood_geometry.getCenterCellStoredMoments())[0].volume() /
          middle_cell.calculateVolume(),
      &initial_reconstruction);

  final_reconstruction =
      reconstructionWithR2P2D(neighborhood_geometry, initial_reconstruction);
  rotateNormalOntoNormal(final_reconstruction[0].normal(),
                         correct_plane_0.normal(), &rotation_angle,
                         &rotation_axis);

  EXPECT_EQ(final_reconstruction.isFlipped(),
            correct_reconstruction.isFlipped());

  EXPECT_EQ(final_reconstruction.getNumberOfPlanes(), 2);

  EXPECT_NEAR(getVolumeFraction(middle_cell, final_reconstruction),
              (neighborhood_geometry.getCenterCellStoredMoments())[0].volume() /
                  middle_cell.calculateVolume(),
              global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);

  std::cout
      << "\n\n"
         "******** Correctness is not automatically checked here ********\n"
      << "Target normal[0]: " << correct_reconstruction[0].normal()[0] << " "
      << correct_reconstruction[0].normal()[1] << " "
      << correct_reconstruction[0].normal()[2] << '\n'
      << "Target normal[1]: " << correct_reconstruction[1].normal()[0] << " "
      << correct_reconstruction[1].normal()[1] << " "
      << correct_reconstruction[1].normal()[2] << '\n'
      << "Target distances: " << correct_reconstruction[0].distance() << " "
      << correct_reconstruction[1].distance() << '\n'
      << "Initial normal[0]: " << initial_reconstruction[0].normal()[0] << " "
      << initial_reconstruction[0].normal()[1] << " "
      << initial_reconstruction[0].normal()[2] << '\n'
      << "Initial normal[1]: " << initial_reconstruction[1].normal()[0] << " "
      << initial_reconstruction[1].normal()[1] << " "
      << initial_reconstruction[1].normal()[2] << '\n'
      << "Initial distances: " << initial_reconstruction[0].distance() << " "
      << initial_reconstruction[1].distance() << '\n'
      << "Found normal[0]:  " << final_reconstruction[0].normal()[0] << " "
      << final_reconstruction[0].normal()[1] << " "
      << final_reconstruction[0].normal()[2] << '\n'
      << "Found normal[1]:  " << final_reconstruction[1].normal()[0] << " "
      << final_reconstruction[1].normal()[1] << " "
      << final_reconstruction[1].normal()[2] << '\n'
      << "Found distances: " << final_reconstruction[0].distance() << " "
      << final_reconstruction[1].distance() << '\n'

      << "The angle between the two shared normals is: "
      << rad2Deg(rotation_angle) << " degrees \n"
      << "******** Check this ^^^^^^^^ for correctness ********\n"
      << "To print out path of optimization (in rotated frame)\n"
      << "change reconstructionWithR2P2D to reconstructionWithR2P2DDebug\n\n";
}

TEST(ReconstructionInterface, R2P_3D2P) {
  Plane correct_plane_0(Normal::normalized(-1.0, 1.0, 1.0), 0.0);
  Plane correct_plane_1(Normal::normalized(1.0, 1.0, 1.0), 0.0);
  PlanarSeparator correct_reconstruction =
      PlanarSeparator::fromTwoPlanes(correct_plane_0, correct_plane_1, 1.0);
  R2PNeighborhood<RectangularCuboid> neighborhood_geometry;
  neighborhood_geometry.resize(27);
  RectangularCuboid cells[27];
  SeparatedMoments<VolumeMoments> svm[27];
  for (int k = -1; k < 2; ++k) {
    for (int j = -1; j < 2; ++j) {
      for (int i = -1; i < 2; ++i) {
        UnsignedIndex_t ind =
            static_cast<UnsignedIndex_t>((k + 1) * 9 + (j + 1) * 3 + i + 1);
        cells[ind] = unit_cell;
        cells[ind].shift(static_cast<double>(i), static_cast<double>(j),
                         static_cast<double>(k));
        svm[ind] = getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
            cells[ind], correct_reconstruction);
        neighborhood_geometry.setMember(ind, &cells[ind], &svm[ind]);
      }
    }
  }
  neighborhood_geometry.setSurfaceArea(
      getReconstructionSurfaceArea(unit_cell, correct_reconstruction));
  neighborhood_geometry.setCenterOfStencil(13);

  // Perturb initial reconstruction then see if it is recovered.
  Normal rotation_axis;
  double rotation_angle;
  Normal shared_normal =
      getSharedNormal(correct_plane_1.normal(), correct_plane_0.normal(),
                      &rotation_angle, &rotation_axis);
  ReferenceFrame ref_frame(rotation_axis,
                           crossProductNormalized(rotation_axis, shared_normal),
                           shared_normal);
  UnitQuaternion perturb_rotate = (UnitQuaternion(deg2Rad(15.0), ref_frame[2]) *
                                   UnitQuaternion(deg2Rad(20.0), ref_frame[1]) *
                                   UnitQuaternion(deg2Rad(18.0), ref_frame[0]));
  Plane initial_plane_0(perturb_rotate * correct_plane_0.normal(), 0.0);
  Plane initial_plane_1(perturb_rotate * correct_plane_1.normal(), 0.0);
  PlanarSeparator initial_reconstruction =
      PlanarSeparator::fromTwoPlanes(initial_plane_0, initial_plane_1, 1.0);
  setDistanceToMatchVolumeFractionPartialFill(
      unit_cell,
      (neighborhood_geometry.getCenterCellStoredMoments())[0].volume(),
      &initial_reconstruction);
  PlanarSeparator final_reconstruction;
  final_reconstruction =
      reconstructionWithR2P3D(neighborhood_geometry, initial_reconstruction);
  rotateNormalOntoNormal(final_reconstruction[0].normal(),
                         correct_plane_0.normal(), &rotation_angle,
                         &rotation_axis);

  EXPECT_EQ(final_reconstruction.isFlipped(),
            correct_reconstruction.isFlipped());

  EXPECT_NEAR(getVolumeFraction(unit_cell, final_reconstruction),
              (neighborhood_geometry.getCenterCellStoredMoments())[0].volume(),
              global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);

  std::cout
      << "\n\n"
         "******** Correctness is not automatically checked here ********\n"
      << "Target normal[0]: " << correct_reconstruction[0].normal()[0] << " "
      << correct_reconstruction[0].normal()[1] << " "
      << correct_reconstruction[0].normal()[2] << '\n'
      << "Target normal[1]: " << correct_reconstruction[1].normal()[0] << " "
      << correct_reconstruction[1].normal()[1] << " "
      << correct_reconstruction[1].normal()[2] << '\n'
      << "Target distances: " << correct_reconstruction[0].distance() << " "
      << correct_reconstruction[1].distance() << '\n'
      << "Found normal[0]:  " << final_reconstruction[0].normal()[0] << " "
      << final_reconstruction[0].normal()[1] << " "
      << final_reconstruction[0].normal()[2] << '\n'
      << "Found normal[1]:  " << final_reconstruction[1].normal()[0] << " "
      << final_reconstruction[1].normal()[1] << " "
      << final_reconstruction[1].normal()[2] << '\n'
      << "Found distances: " << final_reconstruction[0].distance() << " "
      << final_reconstruction[1].distance() << '\n'

      << "The angle between the two shared normals is: "
      << rad2Deg(rotation_angle) << " degrees \n"
      << "******** Check this ^^^^^^^^ for correctness ********\n"
      << "To print out path of optimization (in rotated frame)\n"
      << "change reconstructionWithR2P3D to reconstructionWithR2P3DDebug\n\n";
}

TEST(ReconstructionInterface, ELVIRA_2D) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(
      global_constants::VF_LOW + DBL_EPSILON,
      global_constants::VF_HIGH - DBL_EPSILON);
  RectangularCuboid stencil_cells[9];
  double cellVF[9];
  ELVIRANeighborhood neighborhood_VF;
  neighborhood_VF.resize(9);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      stencil_cells[i + j * 3] = unit_cell;
      stencil_cells[i + j * 3].shift(static_cast<double>(i - 1),
                                     static_cast<double>(j - 1), 0.0);
      neighborhood_VF.setMember(&stencil_cells[i + j * 3], &cellVF[i + j * 3],
                                i - 1, j - 1);
    }
  }
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal correct_normal =
        Normal::normalized(random_normal(eng), random_normal(eng), 0.0);
    double set_VF = random_VF(eng);
    PlanarSeparator correct_reconstruction = PlanarSeparator::fromOnePlane(
        Plane(correct_normal,
              findDistanceOnePlane(unit_cell, set_VF, correct_normal)));
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        cellVF[i + j * 3] =
            getVolumeFraction(stencil_cells[i + j * 3], correct_reconstruction);
      }
    }
    auto found_reconstruction = reconstructionWithELVIRA2D(neighborhood_VF);
    EXPECT_NEAR(found_reconstruction[0].normal()[0],
                correct_reconstruction[0].normal()[0], 1.0e-15);
    EXPECT_NEAR(found_reconstruction[0].normal()[1],
                correct_reconstruction[0].normal()[1], 1.0e-15);
    EXPECT_NEAR(found_reconstruction[0].normal()[2],
                correct_reconstruction[0].normal()[2], 1.0e-15);
    double found_VF = getVolumeFraction(stencil_cells[4], found_reconstruction);
    EXPECT_NEAR(set_VF, found_VF, 1.0e-14);
  }
}

TEST(ReconstructionInterface, ELVIRA_3D) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(
      global_constants::VF_LOW + DBL_EPSILON,
      global_constants::VF_HIGH - DBL_EPSILON);
  RectangularCuboid stencil_cells[27];
  double cellVF[27];
  ELVIRANeighborhood neighborhood_VF;
  neighborhood_VF.resize(27);
  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        stencil_cells[i + j * 3 + k * 9] = unit_cell;
        stencil_cells[i + j * 3 + k * 9].shift(static_cast<double>(i - 1),
                                               static_cast<double>(j - 1),
                                               static_cast<double>(k - 1));
        neighborhood_VF.setMember(&stencil_cells[i + j * 3 + k * 9],
                                  &cellVF[i + j * 3 + k * 9], i - 1, j - 1,
                                  k - 1);
      }
    }
  }
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal correct_normal = Normal::normalized(
        random_normal(eng), random_normal(eng), random_normal(eng));
    double set_VF = random_VF(eng);
    PlanarSeparator correct_reconstruction = PlanarSeparator::fromOnePlane(
        Plane(correct_normal,
              findDistanceOnePlane(unit_cell, set_VF, correct_normal)));

    for (int k = 0; k < 3; ++k) {
      for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
          cellVF[i + j * 3 + k * 9] = getVolumeFraction(
              stencil_cells[i + j * 3 + k * 9], correct_reconstruction);
        }
      }
    }
    auto found_reconstruction = reconstructionWithELVIRA3D(neighborhood_VF);
    EXPECT_NEAR(found_reconstruction[0].normal()[0],
                correct_reconstruction[0].normal()[0], 1.0e-3);
    EXPECT_NEAR(found_reconstruction[0].normal()[1],
                correct_reconstruction[0].normal()[1], 1.0e-3);
    EXPECT_NEAR(found_reconstruction[0].normal()[2],
                correct_reconstruction[0].normal()[2], 1.0e-3);
    double found_VF =
        getVolumeFraction(stencil_cells[13], found_reconstruction);
    EXPECT_NEAR(set_VF, found_VF, 1.0e-14);
  }
}

TEST(ReconstructionInterface, LVIRA_2D) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(
      global_constants::VF_LOW + DBL_EPSILON,
      global_constants::VF_HIGH - DBL_EPSILON);
  RectangularCuboid stencil_cells[9];
  double cellVF[9];
  LVIRANeighborhood<RectangularCuboid> neighborhood_VF;
  neighborhood_VF.resize(9);
  neighborhood_VF.setCenterOfStencil(4);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      stencil_cells[i + j * 3] = unit_cell;
      stencil_cells[i + j * 3].shift(static_cast<double>(i - 1),
                                     static_cast<double>(j - 1), 0.0);
      neighborhood_VF.setMember(static_cast<UnsignedIndex_t>(i + j * 3),
                                &stencil_cells[i + j * 3], &cellVF[i + j * 3]);
    }
  }
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal correct_normal =
        Normal::normalized(random_normal(eng), random_normal(eng), 0.0);
    double set_VF = random_VF(eng);
    PlanarSeparator correct_reconstruction = PlanarSeparator::fromOnePlane(
        Plane(correct_normal,
              findDistanceOnePlane(unit_cell, set_VF, correct_normal)));
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        cellVF[i + j * 3] =
            getVolumeFraction(stencil_cells[i + j * 3], correct_reconstruction);
      }
    }

    auto separated_moments =
        getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
            neighborhood_VF.getCenterCell(), correct_reconstruction);
    auto bary_normal = Normal::fromPtNormalized(
        separated_moments[1].centroid() - separated_moments[0].centroid());
    const double plane_distance =
        bary_normal * neighborhood_VF.getCenterCell().calculateCentroid();
    auto found_reconstruction =
        PlanarSeparator::fromOnePlane(Plane(bary_normal, plane_distance));
    setDistanceToMatchVolumeFractionPartialFill(
        neighborhood_VF.getCenterCell(),
        neighborhood_VF.getCenterCellStoredMoments(), &found_reconstruction);

    found_reconstruction =
        reconstructionWithLVIRA2D(neighborhood_VF, found_reconstruction);
    EXPECT_NEAR(found_reconstruction[0].normal()[0],
                correct_reconstruction[0].normal()[0], 1.0e-3);
    EXPECT_NEAR(found_reconstruction[0].normal()[1],
                correct_reconstruction[0].normal()[1], 1.0e-3);
    EXPECT_NEAR(found_reconstruction[0].normal()[2],
                correct_reconstruction[0].normal()[2], 1.0e-3);
    double found_VF = getVolumeFraction(stencil_cells[4], found_reconstruction);
    EXPECT_NEAR(set_VF, found_VF, 1.0e-14);
  }
}

TEST(ReconstructionInterface, LVIRA_3D) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(
      global_constants::VF_LOW + DBL_EPSILON,
      global_constants::VF_HIGH - DBL_EPSILON);
  RectangularCuboid stencil_cells[27];
  double cellVF[27];
  LVIRANeighborhood<RectangularCuboid> neighborhood_VF;
  neighborhood_VF.resize(27);
  neighborhood_VF.setCenterOfStencil(13);
  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        stencil_cells[i + j * 3 + k * 9] = unit_cell;
        stencil_cells[i + j * 3 + k * 9].shift(static_cast<double>(i - 1),
                                               static_cast<double>(j - 1),
                                               static_cast<double>(k - 1));
        neighborhood_VF.setMember(
            static_cast<UnsignedIndex_t>(i + j * 3 + k * 9),
            &stencil_cells[i + j * 3 + k * 9], &cellVF[i + j * 3 + k * 9]);
      }
    }
  }
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal correct_normal = Normal::normalized(
        random_normal(eng), random_normal(eng), random_normal(eng));
    double set_VF = random_VF(eng);
    PlanarSeparator correct_reconstruction = PlanarSeparator::fromOnePlane(
        Plane(correct_normal,
              findDistanceOnePlane(unit_cell, set_VF, correct_normal)));

    for (int k = 0; k < 3; ++k) {
      for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
          cellVF[i + j * 3 + k * 9] = getVolumeFraction(
              stencil_cells[i + j * 3 + k * 9], correct_reconstruction);
        }
      }
    }
    auto separated_moments =
        getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
            neighborhood_VF.getCenterCell(), correct_reconstruction);
    auto bary_normal = Normal::fromPtNormalized(
        separated_moments[1].centroid() - separated_moments[0].centroid());
    const double plane_distance =
        bary_normal * neighborhood_VF.getCenterCell().calculateCentroid();
    auto found_reconstruction =
        PlanarSeparator::fromOnePlane(Plane(bary_normal, plane_distance));
    setDistanceToMatchVolumeFractionPartialFill(
        neighborhood_VF.getCenterCell(),
        neighborhood_VF.getCenterCellStoredMoments(), &found_reconstruction);

    found_reconstruction =
        reconstructionWithLVIRA3D(neighborhood_VF, found_reconstruction);
    EXPECT_NEAR(found_reconstruction[0].normal()[0],
                correct_reconstruction[0].normal()[0], 1.0e-3);
    EXPECT_NEAR(found_reconstruction[0].normal()[1],
                correct_reconstruction[0].normal()[1], 1.0e-3);
    EXPECT_NEAR(found_reconstruction[0].normal()[2],
                correct_reconstruction[0].normal()[2], 1.0e-3);
    double found_VF =
        getVolumeFraction(stencil_cells[13], found_reconstruction);
    EXPECT_NEAR(set_VF, found_VF, 1.0e-14);
  }
}

TEST(ReconstructionInterface, MOF_2D) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(
      global_constants::VF_LOW + DBL_EPSILON,
      global_constants::VF_HIGH - DBL_EPSILON);
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal correct_normal =
        Normal::normalized(random_normal(eng), random_normal(eng), 0.0);
    double set_VF = random_VF(eng);
    PlanarSeparator correct_reconstruction = PlanarSeparator::fromOnePlane(
        Plane(correct_normal,
              findDistanceOnePlane(unit_cell, set_VF, correct_normal)));
    auto correct_svm =
        getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
            unit_cell, correct_reconstruction);

    auto found_reconstruction = reconstructionWithMOF2D(unit_cell, correct_svm);
    auto found_volume_moments =
        getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
            unit_cell, found_reconstruction);
    EXPECT_NEAR(squaredMagnitude(found_volume_moments[0].centroid() -
                                 correct_svm[0].centroid()),
                0.0, 1.0e-4)
        << "Internal volume: " << set_VF << '\n';
    EXPECT_NEAR(squaredMagnitude(found_volume_moments[1].centroid() -
                                 correct_svm[1].centroid()),
                0.0, 1.0e-4)
        << "External volume: " << 1.0 - set_VF << '\n';
    EXPECT_NEAR(
        std::fabs(found_volume_moments[0].volume() - correct_svm[0].volume()),
        0.0, 1.0e-14);
    EXPECT_NEAR(
        std::fabs(found_volume_moments[1].volume() - correct_svm[1].volume()),
        0.0, 1.0e-14);
  }
}

TEST(ReconstructionInterface, MOF_2D_Triangle) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(
      global_constants::VF_LOW + DBL_EPSILON,
      global_constants::VF_HIGH - DBL_EPSILON);
  Tri base_tri({Pt(0.0, 0.0, 0.0), Pt(1.0, 0.0, 0.0), Pt(0.0, 1.0, 0.0)});
  base_tri.calculateAndSetPlaneOfExistence();
  auto base_tri_volume = base_tri.calculateVolume();
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal correct_normal =
        Normal::normalized(random_normal(eng), random_normal(eng), 0.0);
    double set_VF = random_VF(eng);

    PlanarSeparator correct_reconstruction =
        PlanarSeparator::fromOnePlane(Plane(correct_normal, 0.0));

    setDistanceToMatchVolumeFractionPartialFill(base_tri, set_VF,
                                                &correct_reconstruction);
    auto correct_svm =
        getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
            base_tri, correct_reconstruction);

    auto found_reconstruction = reconstructionWithMOF3D(base_tri, correct_svm);
    auto found_volume_moments =
        getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
            base_tri, found_reconstruction);

    EXPECT_NEAR(squaredMagnitude(found_volume_moments[0].centroid() -
                                 correct_svm[0].centroid()) /
                    std::pow(base_tri_volume, 1.0 / 3.0),
                0.0, 1.0e-4)
        << "Internal volume: " << set_VF << '\n';
    EXPECT_NEAR(squaredMagnitude(found_volume_moments[1].centroid() -
                                 correct_svm[1].centroid()) /
                    std::pow(base_tri_volume, 1.0 / 3.0),
                0.0, 1.0e-4)
        << "External volume: " << 1.0 - set_VF << '\n';
    EXPECT_NEAR(std::fabs((found_volume_moments[0].volume() -
                           correct_svm[0].volume())) /
                    base_tri_volume,
                0.0,
                global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);
    EXPECT_NEAR(std::fabs((found_volume_moments[1].volume() -
                           correct_svm[1].volume())) /
                    base_tri_volume,
                0.0,
                global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);
  }
}

// Generates error map for a random MOF_3D configuration.
// Writes scalarError as a function of two rotation angles.
// Written to file ErrorMap in pm3d format for GNUplot's splot.
// TEST(ReconstructionInterface, MOF_3D_map) {
//  std::random_device
//      rd;  // Get a random seed from the OS entropy device, or whatever
//  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
//                              // generator and seed it with entropy.
//  static const int ncycles = 1;
//  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
//  std::uniform_real_distribution<double> random_VF(
//      global_constants::VF_LOW + DBL_EPSILON,
//      global_constants::VF_HIGH - DBL_EPSILON);
//  for (int cycle = 0; cycle < ncycles; ++cycle) {
//    Normal correct_normal = Normal::normalized(
//        random_normal(eng), random_normal(eng), random_normal(eng));
//    double set_VF = random_VF(eng);
//    std::cout << "Normal " << correct_normal << std::endl;
//    std::cout << "VOF " << set_VF << std::endl;
//    PlanarSeparator correct_reconstruction = PlanarSeparator::fromOnePlane(
//        Plane(correct_normal,
//              findDistanceOnePlane(unit_cell, set_VF, correct_normal)));
//    auto correct_svm =
//        getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
//            unit_cell, correct_reconstruction);
//
//    MOF_3D<RectangularCuboid> mof_solver;
//    mof_solver.setup(
//        CellGroupedMoments<RectangularCuboid,
//        SeparatedMoments<VolumeMoments>>(
//            &unit_cell, &correct_svm),
//        0.5, 0.5);
//    std::ofstream error_map_file("ErrorMap");
//    for (UnsignedIndex_t t2 = 0; t2 < 360; ++t2) {
//      double angle_t2 = deg2Rad(static_cast<double>(t2));
//      for (UnsignedIndex_t t1 = 0; t1 < 360; ++t1) {
//        double angle_t1 = deg2Rad(static_cast<double>(t1));
//        Eigen::Matrix<double, 2, 1> angle_rotation;
//        angle_rotation(0) = angle_t1;
//        angle_rotation(1) = angle_t2;
//        mof_solver.updateGuess(&angle_rotation);
//        double error = mof_solver.calculateScalarError();
//
//        error_map_file << std::scientific << std::setprecision(15) << t1 << "
//        "
//                       << t2 << " " << error << '\n';
//      }
//      error_map_file << '\n';
//    }
//
//    auto found_reconstruction = reconstructionWithMOF3D(unit_cell,
//    correct_svm); auto found_volume_moments =
//        getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
//            unit_cell, found_reconstruction);
//  }
//}

TEST(ReconstructionInterface, MOF_3D) {
  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(
      global_constants::VF_LOW + DBL_EPSILON,
      global_constants::VF_HIGH - DBL_EPSILON);
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Normal correct_normal = Normal::normalized(
        random_normal(eng), random_normal(eng), random_normal(eng));
    double set_VF = random_VF(eng);
    PlanarSeparator correct_reconstruction = PlanarSeparator::fromOnePlane(
        Plane(correct_normal,
              findDistanceOnePlane(unit_cell, set_VF, correct_normal)));
    auto correct_svm =
        getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
            unit_cell, correct_reconstruction);

    auto found_reconstruction = reconstructionWithMOF3D(unit_cell, correct_svm);
    auto found_volume_moments =
        getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
            unit_cell, found_reconstruction);
    EXPECT_NEAR(squaredMagnitude(found_volume_moments[0].centroid() -
                                 correct_svm[0].centroid()),
                0.0, 1.0e-4)
        << "Internal volume: " << set_VF << '\n';
    EXPECT_NEAR(squaredMagnitude(found_volume_moments[1].centroid() -
                                 correct_svm[1].centroid()),
                0.0, 1.0e-4)
        << "External volume: " << 1.0 - set_VF << '\n';
    EXPECT_NEAR(
        std::fabs(found_volume_moments[0].volume() - correct_svm[0].volume()),
        0.0, 1.0e-14);
    EXPECT_NEAR(
        std::fabs(found_volume_moments[1].volume() - correct_svm[1].volume()),
        0.0, 1.0e-14);
  }
}

TEST(ReconstructionInterface, AdvectedPlaneReconstruction) {
  // Construct list that has been previously tested in
  // TEST(Partitioning, PartitionByNormalKMeansHarder)
  ListedVolumeMoments<VolumeMomentsAndNormal> list;
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(-0.2, 0.18, 0.0)),
                                 Normal(-0.2, 1.0, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(0.2, -0.18, 0.0)),
                                 Normal(0.2, -1.0, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(0.2, 0.18, 0.0)),
                                 Normal(0.2, 1.0, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.5, Pt(-0.2, -0.18, 0.0)),
                                 Normal(-0.2, -1.0, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.1, Pt(0.0, -0.25, 0.0)),
                                 Normal(0.0, -1.0, 0.0));
  list += VolumeMomentsAndNormal(VolumeMoments(0.1, Pt(0.0, 0.25, 0.0)),
                                 Normal(0.0, 1.0, 0.0));

  // Normalize all normals.
  for (auto& element : list) {
    element.normal().normalize();
  }
  list.multiplyByVolume();

  // Set expected PlanarSeparator to create neighborhood moments from.
  PlanarSeparator correct_separator = PlanarSeparator::fromTwoPlanes(
      Plane(Normal(0.0, 1.0, 0.0), -0.25), Plane(Normal(0.0, -1.0, 0.0), -0.25),
      -1.0);
  R2PNeighborhood<RectangularCuboid> neighborhood;
  std::array<RectangularCuboid, 27> cells;
  std::array<SeparatedMoments<VolumeMoments>, 27> moments;
  for (UnsignedIndex_t i = 0; i < 3; ++i) {
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      for (UnsignedIndex_t k = 0; k < 3; ++k) {
        cells[i * 9 + j * 3 + k] = unit_cell;
        cells[i * 9 + j * 3 + k].shift(static_cast<double>(i) - 1.0,
                                       static_cast<double>(j) - 1.0,
                                       static_cast<double>(k) - 1.0);
        moments[i * 9 + j * 3 + k] =
            getNormalizedVolumeMoments<SeparatedMoments<VolumeMoments>>(
                cells[i * 9 + j * 3 + k], correct_separator);
        neighborhood.addMember(&cells[i * 9 + j * 3 + k],
                               &moments[i * 9 + j * 3 + k]);
      }
    }
  }
  neighborhood.setCenterOfStencil(13);
  auto computed_separator =
      reconstructionWithAdvectedNormals(list, neighborhood);
  ASSERT_EQ(computed_separator.getNumberOfPlanes(), 2);
  EXPECT_NEAR(computed_separator[0].normal()[0], 0.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[0].normal()[1], 1.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[0].normal()[2], 0.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[0].distance(), -0.25, 1.0e-8);
  EXPECT_NEAR(computed_separator[1].normal()[0], 0.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[1].normal()[1], -1.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[1].normal()[2], 0.0, 1.0e-15);
  EXPECT_NEAR(computed_separator[1].distance(), -0.25, 1.0e-8);
  EXPECT_EQ(computed_separator.isFlipped(), true);
}

}  // namespace
