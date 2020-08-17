// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/interface_reconstruction_methods/plane_distance.h"

#include <algorithm>
#include <cmath>
#include <float.h>
#include <random>

#include "gtest/gtest.h"

#include "src/geometry/general/plane.h"
#include "src/geometry/general/rotations.h"
#include "src/geometry/polyhedrons/hexahedron.h"
#include "src/interface_reconstruction_methods/volume_fraction_matching.h"
#include "src/parameters/constants.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace {

using namespace IRL;

TEST(PlaneDistance, findDistanceOnePlane) {
  std::random_device
      rd; // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd()); // Use the 64-bit Mersenne Twister 19937
                             // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(
      global_constants::VF_LOW + DBL_EPSILON,
      global_constants::VF_HIGH - DBL_EPSILON);
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    double nx = random_normal(eng);
    double ny = random_normal(eng);
    double nz = random_normal(eng);
    double dist = 0.0;
    double VF = random_VF(eng);
    Plane plane(Normal::normalized(nx, ny, nz), dist);
    PlanarSeparator reconstruction = PlanarSeparator::fromOnePlane(plane);
    reconstruction[0].distance() =
        findDistanceOnePlane(unit_cell, VF, reconstruction[0].normal());
    double new_cut_volume = getVolumeFraction(unit_cell, reconstruction);
    EXPECT_NEAR(new_cut_volume, VF, 1.0e-14)
        << "Normal of : " << plane.normal()[0] << " " << plane.normal()[1]
        << " " << plane.normal()[2] << "  and distance of " << dist
        << " which led to a volume of " << VF << '\n';
  }
}

TEST(PlaneDistance, findDistanceOnePlaneHex) {
  std::random_device
      rd; // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd()); // Use the 64-bit Mersenne Twister 19937
                             // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(
      global_constants::VF_LOW + DBL_EPSILON,
      global_constants::VF_HIGH - DBL_EPSILON);
  // Hexahedron that really is a unit cube
  Hexahedron hex({Pt(0.5, -0.5, -0.5), Pt(0.5, 0.5, -0.5), Pt(0.5, 0.5, 0.5),
                  Pt(0.5, -0.5, 0.5), Pt(-0.5, -0.5, -0.5), Pt(-0.5, 0.5, -0.5),
                  Pt(-0.5, 0.5, 0.5), Pt(-0.5, -0.5, 0.5)});

  for (int cycle = 0; cycle < ncycles; ++cycle) {
    double nx = random_normal(eng);
    double ny = random_normal(eng);
    double nz = random_normal(eng);
    double dist = 0.0;
    double VF = random_VF(eng);
    Plane plane(Normal::normalized(nx, ny, nz), dist);
    PlanarSeparator reconstruction = PlanarSeparator::fromOnePlane(plane);
    reconstruction[0].distance() =
        reconstruction[0].normal() * hex.calculateCentroid();
    setDistanceToMatchVolumeFraction(hex, VF, &reconstruction, 1.0e-14);
    double new_cut_volume = getVolumeFraction(hex, reconstruction);
    EXPECT_NEAR(new_cut_volume, VF, 1.0e-14)
        << "Normal of : " << plane.normal()[0] << " " << plane.normal()[1]
        << " " << plane.normal()[2] << "  and distance of " << dist
        << " which led to a volume of " << VF << '\n';
  }
}

TEST(PlaneDistance, findDistanceOnePlaneTet) {
  std::random_device
      rd; // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd()); // Use the 64-bit Mersenne Twister 19937
                             // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_VF(
      global_constants::VF_LOW + DBL_EPSILON,
      global_constants::VF_HIGH - DBL_EPSILON);
  Tet tet({Pt(1.0, 0.0, -0.5), Pt(1.0, 1.0, 0.0), Pt(1.0, 0.0, 0.5),
           Pt(0.0, 0.0, 0.0)});
  for (auto &vertex : tet) {
    vertex += Pt(-10.0, 5.0, 1.5);
  }
  const auto tet_volume = tet.calculateVolume();

  for (int cycle = 0; cycle < ncycles; ++cycle) {
    double nx = random_normal(eng);
    double ny = random_normal(eng);
    double nz = random_normal(eng);
    double dist = 0.0;
    double VF = random_VF(eng);
    Plane plane(Normal::normalized(nx, ny, nz), dist);
    PlanarSeparator reconstruction = PlanarSeparator::fromOnePlane(plane);
    reconstruction[0].distance() =
        findDistanceOnePlane(tet, VF, reconstruction[0].normal());
    double new_cut_volume = getVolumeFraction(tet, reconstruction);
    EXPECT_NEAR(new_cut_volume, VF, 1.0e-14)
        << "Normal of : " << plane.normal()[0] << " " << plane.normal()[1]
        << " " << plane.normal()[2] << "  and distance of " << dist
        << " which led to a volume of " << VF << '\n';
  }
}

TEST(PlaneDistance, findDistanceTwoPlane) {
  std::random_device
      rd; // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd()); // Use the 64-bit Mersenne Twister 19937
                             // generator and seed it with entropy.
  static const int ncycles = 500;
  std::uniform_real_distribution<double> random_normal(-1.0, 1.0);
  std::uniform_real_distribution<double> random_beta(0.0, 2.0 * M_PI);
  std::uniform_real_distribution<double> random_distance(-0.5, 0.5);
  std::uniform_real_distribution<double> random_distance_separation(1.0e-6,
                                                                    0.20);
  std::uniform_real_distribution<double> center_x_location(-100.0, 100.0);
  std::uniform_real_distribution<double> center_y_location(-100.0, 100.0);
  std::uniform_real_distribution<double> center_z_location(-100.0, 100.0);
  std::uniform_real_distribution<double> cell_size_x(1.0e-8, 100.0);
  std::uniform_real_distribution<double> cell_size_y(1.0e-8, 100.0);
  std::uniform_real_distribution<double> cell_size_z(1.0e-8, 100.0);

  std::uniform_real_distribution<double> random_distance_shift(-0.15, 0.15);
  std::uniform_real_distribution<double> random_tolerance(1.0e-14, 1.0e-10);
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    double nx = random_normal(eng);
    double ny = random_normal(eng);
    double nz = random_normal(eng);
    double beta = random_beta(eng);

    double cxl = center_x_location(eng);
    double cyl = center_y_location(eng);
    double czl = center_z_location(eng);
    double csx = cell_size_x(eng);
    double csy = cell_size_y(eng);
    double csz = cell_size_z(eng);
    RectangularCuboid cell = RectangularCuboid::fromBoundingPts(
        Pt(cxl - 0.5 * csx, cyl - 0.5 * csy, czl - 0.5 * csz),
        Pt(cxl + 0.5 * csx, cyl + 0.5 * csy, czl + 0.5 * csz));

    Normal shared_normal = Normal::normalized(nx, ny, nz);
    ReferenceFrame ref_frame = getOrthonormalSystem(shared_normal);
    double dist = random_distance(eng) * (shared_normal * Pt(csx, csy, csz)) +
                  shared_normal * Pt(cxl, cyl, czl);
    double dist_separation =
        random_distance_separation(eng) * std::pow(csx * csy * csz, 1.0 / 3.0);
    double dist_shift =
        random_distance_shift(eng) * std::pow(csx * csy * csz, 1.0 / 3.0);
    double tolerance = random_tolerance(eng);

    Plane plane_0(UnitQuaternion(beta, ref_frame[0]) * ref_frame[2], dist);
    Plane plane_1(UnitQuaternion(-beta, ref_frame[0]) * ref_frame[2],
                  dist + std::copysign(dist_separation, dist));
    PlanarSeparator reconstruction =
        PlanarSeparator::fromTwoPlanes(plane_0, plane_1, 1.0);
    if (beta > M_PI) {
      reconstruction.flipCutting();
    }
    double cut_volume = getVolumeFraction(cell, reconstruction);
    std::cout << "cut volume " << cut_volume << std::endl;
    SmallVector<double, 2> distances;
    distances.resize(2);
    distances[0] = plane_0.distance() + dist_shift;
    distances[1] = plane_1.distance() - dist_shift;
    reconstruction.setDistances(distances);
    setDistanceToMatchVolumeFraction(cell, cut_volume, &reconstruction,
                                     tolerance);
    double new_cut_volume = getVolumeFraction(cell, reconstruction);
    EXPECT_NEAR(new_cut_volume, cut_volume, tolerance)
        << "Plane 0 : " << plane_0.normal()[0] << " " << plane_0.normal()[1]
        << " " << plane_0.normal()[2] << " " << plane_0.distance() << "\n"
        << "Plane 1 : " << plane_1.normal()[0] << " " << plane_1.normal()[1]
        << " " << plane_1.normal()[2] << " " << plane_1.distance() << "\n"
        << "Given dist : " << distances[0] << " " << distances[1] << '\n'
        << "Cell: \n"
        << cell[4].x() << " " << cell[4].y() << " " << cell[4].z() << '\n'
        << cell[2].x() << " " << cell[2].y() << " " << cell[2].z() << '\n';
  }
}

TEST(PlaneDistance, hangingTwoPlane) {
  RectangularCuboid cell = RectangularCuboid::fromBoundingPts(
      Pt(-5.0e-5, -5.0e-5, -5.0e-5), Pt(5.0e-5, 5.0e-5, 5.0e-5));

  double cell_vol = std::pow(1.0e-4, 3);
  setVolumeFractionBounds(1.0e-8);
  setVolumeFractionTolerance(1.0e-12);
  setMinimumVolumeToTrack(cell_vol * 1.0e-15);
  setMinimumSurfaceAreaToTrack(cell_vol * 1.0e-15);

  double correct_volume_fraction = 5.43141744433826e-07;
  Plane plane_0 =
      Plane(Normal(-0.410773177948303, -0.732471255731464, -0.542909988677152),
            -2.12391961534382e-05);
  Plane plane_1 =
      Plane(Normal(0.408062862804978, 0.616743458406921, 0.673136098060174),
            2.12362107479903e-05);
  PlanarSeparator intial_reconstruction =
      PlanarSeparator::fromTwoPlanes(plane_0, plane_1, 1.0);
  setDistanceToMatchVolumeFractionPartialFill(cell, correct_volume_fraction,
                                              &intial_reconstruction);
}

TEST(PlaneDistance, findDistanceThreePlane) {
  // Simple test for now
  RectangularCuboid cell = unit_cell;
  PlanarSeparator open_box;
  open_box.addPlane(Plane(Normal(0.0, 1.0, 0.0), 0.0));
  open_box.addPlane(Plane(Normal(-1.0, 0.0, 0.0), 0.1));
  open_box.addPlane(Plane(Normal(1.0, 0.0, 0.0), 0.1));
  EXPECT_NEAR(getVolumeMoments<Volume>(cell, open_box), 0.1, 1.0e-15);
  double target_volume = 0.5;
  IterativeSolverForDistance<ReconstructionDefaultCuttingMethod,
                             RectangularCuboid, 3>
      solver(cell, target_volume,
             global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE,
             open_box);
  open_box.setDistances(solver.getDistances());
  double new_cut_volume = getVolumeMoments<Volume>(cell, open_box);
  EXPECT_NEAR(new_cut_volume, target_volume,
              global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);

  PlanarSeparator triangle;
  Normal tmp_normal =
      Normal::normalized(-0.5 * std::sqrt(2.0), 0.5 * std::sqrt(2.0), 0.0);
  triangle.addPlane(Plane(tmp_normal, tmp_normal * Pt(-0.5, 0.0, 0.0)));
  tmp_normal[0] = -tmp_normal[0];
  triangle.addPlane(Plane(tmp_normal, tmp_normal * Pt(0.5, 0.0, 0.0)));
  tmp_normal = Normal(0.0, -1.0, 0.0);
  triangle.addPlane(Plane(tmp_normal, 0.0));
  EXPECT_NEAR(getVolumeMoments<Volume>(cell, triangle), 0.25, 1.0e-15);
  target_volume = 0.75;
  solver = IterativeSolverForDistance<ReconstructionDefaultCuttingMethod,
                                      RectangularCuboid, 3>(
      cell, target_volume,
      global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE, open_box);
  open_box.setDistances(solver.getDistances());
  new_cut_volume = getVolumeMoments<Volume>(cell, open_box);
  EXPECT_NEAR(new_cut_volume, target_volume,
              global_constants::TWO_PLANE_DISTANCE_VOLUME_FRACTION_TOLERANCE);
}

} // namespace
