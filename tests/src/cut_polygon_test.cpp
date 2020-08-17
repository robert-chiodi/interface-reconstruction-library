// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/generic_cutting/cut_polygon.h"

#include <float.h>
#include <algorithm>
#include <cmath>
#include <random>
#include "src/planar_reconstruction/planar_separator.h"

#include "gtest/gtest.h"

#include "src/geometry/general/plane.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/helpers/geometric_cutting_helpers.h"
#include "src/helpers/helper.h"
#include "src/helpers/mymath.h"

namespace {

using namespace IRL;

TEST(CutPolygon, cutPlaneByHexahedron) {
  Plane mid_plane(Normal(0.0, 1.0, 0.0), 0.0);
  Polygon mid_plane_poly = cutPlaneByHexahedron<Polygon>(unit_cell, mid_plane);
  mid_plane_poly.calculateAndSetPlaneOfExistence();
  EXPECT_DOUBLE_EQ(mid_plane_poly.calculateVolume(), 1.0);

  Plane diagonal_plane(Normal::normalized(1.0, 1.0, 0.0), 0.0);
  Polygon diagonal_plane_poly =
      cutPlaneByHexahedron<Polygon>(unit_cell, diagonal_plane);
  diagonal_plane_poly.calculateAndSetPlaneOfExistence();
  EXPECT_DOUBLE_EQ(diagonal_plane_poly.calculateVolume(), std::sqrt(2.0));

  std::random_device
      rd;  // Get a random seed from the OS entropy device, or whatever
  std::mt19937_64 eng(rd());  // Use the 64-bit Mersenne Twister 19937
                              // generator and seed it with entropy.

  // Test rectangles (since will have 4 intersections)
  static const int ncycles = 50;
  std::uniform_real_distribution<double> random_dist(-0.5, 0.5);
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Pt tri1(random_dist(eng), -0.5, -0.5);
    Pt tri2(random_dist(eng), 0.5, -0.5);
    Pt tri3(random_dist(eng), -0.5, 0.5);
    tri2[0] = tri1[0];
    Normal plane_normal = crossProductNormalized(Normal::fromPt(tri2 - tri1),
                                                 Normal::fromPt(tri3 - tri1));
    double dist = plane_normal * tri1;
    Plane cut_plane(plane_normal, dist);
    Polygon polygon = cutPlaneByHexahedron<Polygon>(unit_cell, cut_plane);
    polygon.calculateAndSetPlaneOfExistence();
    EXPECT_DOUBLE_EQ(
        polygon.calculateVolume(),
        distanceBetweenPts(tri1, tri2) * distanceBetweenPts(tri1, tri3));
  }

  // Test Tris (since will have 3 intersections)
  for (int cycle = 0; cycle < ncycles; ++cycle) {
    Pt tri1(-0.5, random_dist(eng), -0.5);
    Pt tri2(-0.5, -0.5, random_dist(eng));
    Pt tri3(random_dist(eng), -0.5, -0.5);
    Normal plane_normal = Normal::fromPtNormalized(
        crossProduct(Pt(tri2 - tri1), Pt(tri3 - tri1)));
    double area =
        0.5 * magnitude(crossProduct(Pt(tri2 - tri1), Pt(tri3 - tri1)));
    double dist = plane_normal * tri1;
    Plane cut_plane(plane_normal, dist);
    Polygon polygon = cutPlaneByHexahedron<Polygon>(unit_cell, cut_plane);
    EXPECT_NEAR(polygon.calculateVolume(), area, 1.0e-13)
        << "Tri Pt1 : " << tri1[0] << " " << tri1[1] << " " << tri1[2] << '\n'
        << "Tri Pt1 : " << tri2[0] << " " << tri2[1] << " " << tri2[2] << '\n'
        << "Tri Pt1 : " << tri3[0] << " " << tri3[1] << " " << tri3[2] << '\n';
  }
}

TEST(CutPolygon, cutPlaneByTet) {
  Tet tet({Pt(1.0, 0.0, -0.5), Pt(1.0, 1.0, 0.0), Pt(1.0, 0.0, 0.5),
           Pt(0.0, 0.0, 0.0)});
  auto tet_volume = tet.calculateVolume();
  ASSERT_TRUE(tet_volume > 0.0);

  std::array<double, 4> distance_to_vertices;

  // Now test all the cutting cases and make sure resulting volume for
  // under/over plane is positive.
  PlanarSeparator planar_reconstruction =
      PlanarSeparator::fromOnePlane(Plane(Normal(0.0, 0.0, 0.0), 0.0));

  double expected_area;
  double area;
  Polygon poly;

  // case 0
  Normal plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 2.0);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 0);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  expected_area = 0.0;
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 1
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 1);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[1], distance_to_vertices[1]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[2], distance_to_vertices[2]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[3], distance_to_vertices[3]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 2
  plane_normal = Normal(0.0, 1.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 2);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[1], distance_to_vertices[1]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[2], distance_to_vertices[2]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 3
  plane_normal = Normal::normalized(1.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(plane_normal, tet[2] * plane_normal * 1.1);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 3);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[2], distance_to_vertices[2]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[2], distance_to_vertices[2],
                                          tet[0], distance_to_vertices[0]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[3], distance_to_vertices[3]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 4
  plane_normal = Normal::normalized(0.0, 0.0, 1.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.1);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 4);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[2], distance_to_vertices[2]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[2], distance_to_vertices[2]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[2], distance_to_vertices[2],
                                          tet[3], distance_to_vertices[3]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 5
  plane_normal = Normal::normalized(1.0, -1.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, tet[2] * plane_normal * 0.9);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 5);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[2], distance_to_vertices[2]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[2], distance_to_vertices[2],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[3], distance_to_vertices[3],
                                          tet[0], distance_to_vertices[0]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[0], distance_to_vertices[0]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 6
  plane_normal = Normal::normalized(1.0, 0.0, 1.0);
  planar_reconstruction[0] = Plane(plane_normal, tet[1] * plane_normal * 0.9);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 6);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[1], distance_to_vertices[1]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[2], distance_to_vertices[2],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[2], distance_to_vertices[2],
                                          tet[0], distance_to_vertices[0]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 7
  plane_normal = Normal::normalized(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(plane_normal, 0.5);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 7);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[2], distance_to_vertices[2],
                                          tet[3], distance_to_vertices[3]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 15
  plane_normal = Normal(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -2.0);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 15);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  expected_area = 0.0;
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 14
  plane_normal = Normal(0.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 14);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[1], distance_to_vertices[1]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[2], distance_to_vertices[2]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 13
  plane_normal = Normal(0.0, 1.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 13);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[1], distance_to_vertices[1]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[2], distance_to_vertices[2]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[3], distance_to_vertices[3]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 12
  plane_normal = Normal::normalized(1.0, 0.0, -1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -tet[2] * plane_normal * 1.1);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 12);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[2], distance_to_vertices[2],
                                          tet[0], distance_to_vertices[0]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[2], distance_to_vertices[2]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 11
  plane_normal = Normal::normalized(0.0, 0.0, 1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.1);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 11);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[2], distance_to_vertices[2]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[2], distance_to_vertices[2],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[2], distance_to_vertices[2]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 10
  plane_normal = Normal::normalized(1.0, -1.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -tet[2] * plane_normal * 0.9);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 10);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[2], distance_to_vertices[2]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[0], distance_to_vertices[0]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[3], distance_to_vertices[3],
                                          tet[0], distance_to_vertices[0]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[2], distance_to_vertices[2],
                                          tet[3], distance_to_vertices[3]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 9
  plane_normal = Normal::normalized(1.0, 0.0, 1.0);
  planar_reconstruction[0] = Plane(-plane_normal, -tet[1] * plane_normal * 0.9);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 9);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[1], distance_to_vertices[1]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[2], distance_to_vertices[2],
                                          tet[0], distance_to_vertices[0]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[2], distance_to_vertices[2],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[3], distance_to_vertices[3]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);

  // case 8
  plane_normal = Normal::normalized(1.0, 0.0, 0.0);
  planar_reconstruction[0] = Plane(-plane_normal, -0.5);
  signedDistanceToVertices(tet, planar_reconstruction[0],
                           &distance_to_vertices);
  ASSERT_EQ(getGeometricCaseId(distance_to_vertices), 8);
  area = getReconstructionSurfaceArea(tet, planar_reconstruction);
  poly.setNumberOfVertices(0);
  poly.setPlaneOfExistence(planar_reconstruction[0]);
  poly.addVertex(Pt::fromEdgeIntersection(tet[0], distance_to_vertices[0],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[2], distance_to_vertices[2],
                                          tet[3], distance_to_vertices[3]));
  poly.addVertex(Pt::fromEdgeIntersection(tet[1], distance_to_vertices[1],
                                          tet[3], distance_to_vertices[3]));
  expected_area = poly.calculateVolume();
  ASSERT_NEAR(area, expected_area, 1.0e-14);
}

TEST(CutPolygon, signedDistanceToPolygonVertices) {
  Polygon rectangle;
  rectangle.addVertex(Pt(0.0, 0.0, 0.0));
  rectangle.addVertex(Pt(1.0, 0.0, 0.0));
  rectangle.addVertex(Pt(1.0, 1.0, 0.0));
  rectangle.addVertex(Pt(0.0, 1.0, 0.0));
  Plane cutting_plane(Normal(0.0, 1.0, 0.0), 0.5);
  std::array<double, 6> distance_to_vertices;
  signedDistanceToPolygonVertices(rectangle, cutting_plane, 1.0,
                                  &distance_to_vertices);
  EXPECT_DOUBLE_EQ(distance_to_vertices[0], -0.5);
  EXPECT_DOUBLE_EQ(distance_to_vertices[1], -0.5);
  EXPECT_DOUBLE_EQ(distance_to_vertices[2], 0.5);
  EXPECT_DOUBLE_EQ(distance_to_vertices[3], 0.5);
  EXPECT_DOUBLE_EQ(distance_to_vertices[4], -0.5);
  EXPECT_DOUBLE_EQ(distance_to_vertices[5], -0.5);

  signedDistanceToPolygonVertices(rectangle, cutting_plane, -1.0,
                                  &distance_to_vertices);
  EXPECT_DOUBLE_EQ(distance_to_vertices[0], 0.5);
  EXPECT_DOUBLE_EQ(distance_to_vertices[1], 0.5);
  EXPECT_DOUBLE_EQ(distance_to_vertices[2], -0.5);
  EXPECT_DOUBLE_EQ(distance_to_vertices[3], -0.5);
  EXPECT_DOUBLE_EQ(distance_to_vertices[4], 0.5);
  EXPECT_DOUBLE_EQ(distance_to_vertices[5], 0.5);

  Polygon Tri;
  Tri.addVertex(Pt(0.0, 0.0, 0.0));
  Tri.addVertex(Pt(1.0, 0.0, 0.0));
  Tri.addVertex(Pt(0.0, 1.0, 0.0));
  signedDistanceToPolygonVertices(Tri, cutting_plane, 1.0,
                                  &distance_to_vertices);
  EXPECT_DOUBLE_EQ(distance_to_vertices[0], -0.5);
  EXPECT_DOUBLE_EQ(distance_to_vertices[1], -0.5);
  EXPECT_DOUBLE_EQ(distance_to_vertices[2], 0.5);
  for (UnsignedIndex_t n = 3; n < 6; ++n) {
    EXPECT_DOUBLE_EQ(distance_to_vertices[0], -0.5);
  }
}

TEST(CutPolygon, cutPolygonByPlaneInPlace) {
  Polygon initial_tri;
  initial_tri.addVertex(Pt(0.0, 0.0, 0.0));
  initial_tri.addVertex(Pt(1.0, 0.0, 0.0));
  initial_tri.addVertex(Pt(0.0, 1.0, 0.0));
  Normal plane_normal(1.0, 0.0, 0.0);
  Plane cutting_plane(plane_normal, plane_normal * Pt(0.75, 0.0, 0.0));
  Polygon cut_tri = initial_tri;
  cutPolygonByPlaneInPlace(&cut_tri, cutting_plane, 1.0);
  ASSERT_EQ(cut_tri.getNumberOfVertices(), 4);
  EXPECT_NEAR(cut_tri[0].x(), 0.0, 1.0e-14);
  EXPECT_NEAR(cut_tri[0].y(), 0.0, 1.0e-14);
  EXPECT_NEAR(cut_tri[0].z(), 0.0, 1.0e-14);
  EXPECT_NEAR(cut_tri[1].x(), 0.75, 1.0e-14);
  EXPECT_NEAR(cut_tri[1].y(), 0.0, 1.0e-14);
  EXPECT_NEAR(cut_tri[1].z(), 0.0, 1.0e-14);
  EXPECT_NEAR(cut_tri[2].x(), 0.75, 1.0e-14);
  EXPECT_NEAR(cut_tri[2].y(), 0.25, 1.0e-14);
  EXPECT_NEAR(cut_tri[2].z(), 0.0, 1.0e-14);
  EXPECT_NEAR(cut_tri[3].x(), 0.0, 1.0e-14);
  EXPECT_NEAR(cut_tri[3].y(), 1.0, 1.0e-14);
  EXPECT_NEAR(cut_tri[3].z(), 0.0, 1.0e-14);
}

TEST(CutPolygon, cutPolygonByPlane) {
  Plane orig_plane(Normal::normalized(-1.0, 1.0, 0.0), 0.0);
  Plane cutting_plane(Normal::normalized(1.0, 1.0, 0.0), 0.0);

  Polygon orig_polygon = cutPlaneByHexahedron<Polygon>(unit_cell, orig_plane);
  orig_polygon = cutPolygonByPlane<Polygon>(orig_polygon, cutting_plane, 1.0);
  EXPECT_DOUBLE_EQ(orig_polygon[0].x(), -0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[0].y(), -0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[0].z(), -0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[1].x(), -0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[1].y(), -0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[1].z(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[2].x(), -0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[2].y(), -0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[2].z(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[3].x(), 0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[3].y(), 0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[3].z(), -0.5);

  orig_polygon = cutPlaneByHexahedron<Polygon>(unit_cell, orig_plane);
  orig_polygon = cutPolygonByPlane<Polygon>(orig_polygon, cutting_plane, -1.0);
  EXPECT_DOUBLE_EQ(orig_polygon[0].x(), 0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[0].y(), 0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[0].z(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[1].x(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[1].y(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[1].z(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[2].x(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[2].y(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[2].z(), -0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[3].x(), -0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[3].y(), -0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[3].z(), -0.5);

  orig_polygon = cutPlaneByHexahedron<Polygon>(unit_cell, orig_plane);
  cutting_plane.distance() = 15.0;
  Polygon same_polygon =
      cutPolygonByPlane<Polygon>(orig_polygon, cutting_plane, 1.0);
  EXPECT_EQ(same_polygon.getNumberOfVertices(), 4);
  for (UnsignedIndex_t n = 0; n < 4; ++n) {
    EXPECT_EQ(same_polygon[n].x(), orig_polygon[n].x());
    EXPECT_EQ(same_polygon[n].y(), orig_polygon[n].y());
    EXPECT_EQ(same_polygon[n].z(), orig_polygon[n].z());
  }

  orig_polygon = cutPlaneByHexahedron<Polygon>(unit_cell, orig_plane);
  cutting_plane.distance() = -15.0;
  Polygon no_polygon =
      cutPolygonByPlane<Polygon>(orig_polygon, cutting_plane, 1.0);
  EXPECT_EQ(no_polygon.getNumberOfVertices(), 0);
}

TEST(CutPolygon, getPlanePolygonFromReconstruction) {
  Plane orig_plane(Normal::normalized(-1.0, 1.0, 0.0), 0.0);
  Plane cutting_plane(Normal::normalized(1.0, 1.0, 0.0), 0.0);
  PlanarSeparator reconstruction =
      PlanarSeparator::fromTwoPlanes(orig_plane, cutting_plane, 1.0);

  Polygon orig_polygon = getPlanePolygonFromReconstruction<Polygon>(
      unit_cell, reconstruction, reconstruction[0]);
  EXPECT_DOUBLE_EQ(orig_polygon[0].x(), -0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[0].y(), -0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[0].z(), -0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[1].x(), -0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[1].y(), -0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[1].z(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[2].x(), -0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[2].y(), -0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[2].z(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[3].x(), 0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[3].y(), 0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[3].z(), -0.5);

  reconstruction.flipCutting();
  orig_polygon = getPlanePolygonFromReconstruction<Polygon>(
      unit_cell, reconstruction, reconstruction[0]);
  EXPECT_DOUBLE_EQ(orig_polygon[0].x(), 0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[0].y(), 0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[0].z(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[1].x(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[1].y(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[1].z(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[2].x(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[2].y(), 0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[2].z(), -0.5);
  EXPECT_DOUBLE_EQ(orig_polygon[3].x(), -0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[3].y(), -0.0);
  EXPECT_DOUBLE_EQ(orig_polygon[3].z(), -0.5);

  orig_polygon = cutPlaneByHexahedron<Polygon>(unit_cell, orig_plane);
  reconstruction.doNotFlipCutting();
  reconstruction[1].distance() = 15.0;
  Polygon same_polygon = getPlanePolygonFromReconstruction<Polygon>(
      unit_cell, reconstruction, reconstruction[0]);
  EXPECT_EQ(same_polygon.getNumberOfVertices(), 4);
  for (UnsignedIndex_t n = 0; n < 4; ++n) {
    EXPECT_EQ(same_polygon[n].x(), orig_polygon[n].x());
    EXPECT_EQ(same_polygon[n].y(), orig_polygon[n].y());
    EXPECT_EQ(same_polygon[n].z(), orig_polygon[n].z());
  }

  orig_polygon = cutPlaneByHexahedron<Polygon>(unit_cell, orig_plane);
  reconstruction[1].distance() = -15.0;
  Polygon no_polygon = getPlanePolygonFromReconstruction<Polygon>(
      unit_cell, reconstruction, reconstruction[0]);
  EXPECT_EQ(no_polygon.getNumberOfVertices(), 0);
}

TEST(CutPolygon, getReconstructionSurfaceArea) {
  Plane plane_0(Normal::normalized(-1.0, 1.0, 0.0), 0.0);
  PlanarSeparator single_recon = PlanarSeparator::fromOnePlane(plane_0);
  EXPECT_DOUBLE_EQ(getReconstructionSurfaceArea(unit_cell, single_recon),
                   std::sqrt(2.0));

  plane_0.distance() = 50.0;
  PlanarSeparator all_liq_recon = PlanarSeparator::fromOnePlane(plane_0);
  EXPECT_DOUBLE_EQ(getReconstructionSurfaceArea(unit_cell, all_liq_recon), 0.0);

  plane_0.distance() = -50.0;
  PlanarSeparator all_gas_recon = PlanarSeparator::fromOnePlane(plane_0);
  EXPECT_DOUBLE_EQ(getReconstructionSurfaceArea(unit_cell, all_gas_recon), 0.0);

  plane_0.normal() = Normal::normalized(-1.0, 1.0, 0.0);
  plane_0.distance() = 0.0;

  Plane plane_1(Normal::normalized(1.0, 1.0, 0.0), 0.0);
  PlanarSeparator Tri_recon =
      PlanarSeparator::fromTwoPlanes(plane_0, plane_1, 1.0);
  EXPECT_DOUBLE_EQ(getReconstructionSurfaceArea(unit_cell, Tri_recon),
                   2.0 * std::sqrt(0.5));

  plane_0.normal() = Normal(0.0, 0.0, 1.0);
  plane_0.distance() = -0.25;
  plane_1.normal() = Normal(0.0, 1.0, 0.0);
  plane_1.distance() = -0.25;
  PlanarSeparator corner_recon =
      PlanarSeparator::fromTwoPlanes(plane_0, plane_1, 1.0);
  EXPECT_DOUBLE_EQ(getReconstructionSurfaceArea(unit_cell, corner_recon), 0.5);

  plane_0.normal() = Normal(0.0, 0.0, -1.0);
  plane_0.distance() = -0.25;
  plane_1.normal() = Normal(0.0, -1.0, 0.0);
  plane_1.distance() = -0.25;
  PlanarSeparator inverted_corner_recon =
      PlanarSeparator::fromTwoPlanes(plane_0, plane_1, -1.0);
  EXPECT_DOUBLE_EQ(
      getReconstructionSurfaceArea(unit_cell, inverted_corner_recon), 1.5);

  plane_0.normal() = Normal(0.0, 0.0, 1.0);
  plane_0.distance() = 0.25;
  plane_1.normal() = Normal(0.0, 0.0, -1.0);
  plane_1.distance() = 0.25;
  PlanarSeparator sheet_recon =
      PlanarSeparator::fromTwoPlanes(plane_0, plane_1, 1.0);
  EXPECT_DOUBLE_EQ(getReconstructionSurfaceArea(unit_cell, sheet_recon), 2.0);

  plane_0.normal() = Normal(0.0, 0.0, -1.0);
  plane_1.distance() = -0.25;
  plane_1.normal() = Normal(0.0, 0.0, 1.0);
  plane_1.distance() = -0.25;
  PlanarSeparator gas_sheet_recon =
      PlanarSeparator::fromTwoPlanes(plane_0, plane_1, -1.0);
  EXPECT_DOUBLE_EQ(getReconstructionSurfaceArea(unit_cell, gas_sheet_recon),
                   2.0);
}

}  // namespace
