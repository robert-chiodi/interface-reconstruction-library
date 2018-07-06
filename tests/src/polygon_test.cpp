// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/polygons/polygon.h"

#include <random>

#include "gtest/gtest.h"

#include "src/generic_cutting/cut_polygon.h"
#include "src/generic_cutting/generic_cutting.h"
#include "src/geometry/general/normal.h"
#include "src/geometry/general/rotations.h"
#include "src/geometry/polygons/tri.h"
#include "src/helpers/helper.h"
#include "src/planar_reconstruction/planar_localizer.h"

namespace {

using namespace IRL;

TEST(Polygon, Polygon) {
  // Check construction and access
  Polygon polygon;
  polygon.setNumberOfVertices(7);
  for (UnsignedIndex_t i = 0; i < polygon.getNumberOfVertices(); ++i) {
    polygon[i] = Pt(i + 1, i + 1 * 8, i + 1 * 2 * 8);
  }
  for (UnsignedIndex_t i = 0; i < 7; ++i) {
    Pt vertex = polygon[i];
    EXPECT_DOUBLE_EQ(vertex.x(), i + 1);
    EXPECT_DOUBLE_EQ(vertex.y(), i + 1 * 8);
    EXPECT_DOUBLE_EQ(vertex.z(), i + 1 * 2 * 8);
  }

  Polygon added_polygon;
  for (UnsignedIndex_t i = 0; i < 7; ++i) {
    added_polygon.addVertex(Pt(i + 1, i + 1 * 8, i + 1 * 2 * 8));
  }

  Polygon short_polygon;
  for (UnsignedIndex_t i = 0; i < 4; ++i) {
    short_polygon.addVertex(Pt(i + 1, i + 1 * 8, i + 1 * 2 * 8));
  }

  // Calculate surface area of polygons
  Polygon rectangle;
  rectangle.addVertex(Pt(0.0, 0.0, 0.0));
  rectangle.addVertex(Pt(2.0, 0.0, 0.0));
  rectangle.addVertex(Pt(2.0, 2.0, 0.0));
  rectangle.addVertex(Pt(0.0, 2.0, 0.0));
  rectangle.calculateAndSetPlaneOfExistence();
  EXPECT_DOUBLE_EQ(rectangle.calculateVolume(), 4.0);

  Polygon triangle;
  triangle.addVertex(Pt(0.0, 1.0, 1.0));
  triangle.addVertex(Pt(0.0, 0.0, 0.0));
  triangle.addVertex(Pt(0.0, 1.0, -1.0));
  triangle.calculateAndSetPlaneOfExistence();
  EXPECT_DOUBLE_EQ(triangle.calculateVolume(), 1.0);

  // Check transfer of plane of existence
  auto plane_of_existence = rectangle.getPlaneOfExistence();
  EXPECT_NEAR(plane_of_existence.normal()[0], 0.0, 1.0e-15);
  EXPECT_NEAR(plane_of_existence.normal()[1], 0.0, 1.0e-15);
  EXPECT_NEAR(plane_of_existence.normal()[2], 1.0, 1.0e-15);
  EXPECT_NEAR(plane_of_existence.distance(), 0.0, 1.0e-15);
  auto divided_rectangle = DividedPolygon::fromPolygon(rectangle);
  plane_of_existence = divided_rectangle.getPlaneOfExistence();
  EXPECT_NEAR(plane_of_existence.normal()[0], 0.0, 1.0e-15);
  EXPECT_NEAR(plane_of_existence.normal()[1], 0.0, 1.0e-15);
  EXPECT_NEAR(plane_of_existence.normal()[2], 1.0, 1.0e-15);
  EXPECT_NEAR(plane_of_existence.distance(), 0.0, 1.0e-15);

  for (UnsignedIndex_t n = 0;
       n < divided_rectangle.getNumberOfSimplicesInDecomposition(); ++n) {
    auto proxy_tri = divided_rectangle.getSimplexFromDecomposition(n);
    plane_of_existence = proxy_tri.getPlaneOfExistence();
    EXPECT_NEAR(plane_of_existence.normal()[0], 0.0, 1.0e-15);
    EXPECT_NEAR(plane_of_existence.normal()[1], 0.0, 1.0e-15);
    EXPECT_NEAR(plane_of_existence.normal()[2], 1.0, 1.0e-15);
    EXPECT_NEAR(plane_of_existence.distance(), 0.0, 1.0e-15);
    Tri stored_tri = static_cast<Tri>(proxy_tri);
    plane_of_existence = stored_tri.getPlaneOfExistence();
    EXPECT_NEAR(plane_of_existence.normal()[0], 0.0, 1.0e-15);
    EXPECT_NEAR(plane_of_existence.normal()[1], 0.0, 1.0e-15);
    EXPECT_NEAR(plane_of_existence.normal()[2], 1.0, 1.0e-15);
    EXPECT_NEAR(plane_of_existence.distance(), 0.0, 1.0e-15);
  }
}

TEST(Polygon, getLocalizer) {
  Polygon smaller_quad;
  smaller_quad.addVertex(Pt(0.25, 0.0, 0.0));
  smaller_quad.addVertex(Pt(1.25, 0.0, 0.0));
  smaller_quad.addVertex(Pt(1.25, 1.0, 0.0));
  smaller_quad.addVertex(Pt(0.25, 1.0, 0.0));
  smaller_quad.calculateAndSetPlaneOfExistence();

  Polygon larger_quad = smaller_quad;
  auto initial_centroid = larger_quad.calculateCentroid();
  for (auto& vertex : larger_quad) {
    vertex += 2.0 * (vertex - initial_centroid);
  }
  larger_quad.calculateAndSetPlaneOfExistence();
  PlanarLocalizer smaller_quad_localizer = smaller_quad.getLocalizer();
  auto localized_volume_moments = getNormalizedVolumeMoments<VolumeMoments>(
      larger_quad, smaller_quad_localizer);
  auto smaller_quad_moments = smaller_quad.calculateMoments();
  smaller_quad_moments.normalizeByVolume();
  EXPECT_NEAR(smaller_quad_moments.volume(), localized_volume_moments.volume(),
              1.0e-14);
  EXPECT_NEAR(smaller_quad_moments.centroid()[0],
              localized_volume_moments.centroid()[0], 1.0e-14);
  EXPECT_NEAR(smaller_quad_moments.centroid()[1],
              localized_volume_moments.centroid()[1], 1.0e-14);
  EXPECT_NEAR(smaller_quad_moments.centroid()[2],
              localized_volume_moments.centroid()[2], 1.0e-14);
}

TEST(Polygon, calculateNearestPtOnSurface) {
  Polygon quad;
  quad.addVertex(Pt(0.25, 0.0, 0.0));
  quad.addVertex(Pt(1.25, 0.0, 0.0));
  quad.addVertex(Pt(1.25, 1.0, 0.0));
  quad.addVertex(Pt(0.25, 1.0, 0.0));
  Pt pt_to_project(0.0, 0.0, 0.0);
  Pt nearest_pt = quad.calculateNearestPtOnSurface(pt_to_project);
  // Should project to pt
  EXPECT_NEAR(nearest_pt[0], quad[0][0], 1.0e-14);
  EXPECT_NEAR(nearest_pt[1], quad[0][1], 1.0e-14);
  EXPECT_NEAR(nearest_pt[2], quad[0][2], 1.0e-14);

  // Should project to edge
  pt_to_project = Pt(0.75, 2.0, 0.0);
  nearest_pt = quad.calculateNearestPtOnSurface(pt_to_project);
  EXPECT_NEAR(nearest_pt[0], 0.75, 1.0e-14);
  EXPECT_NEAR(nearest_pt[1], 1.0, 1.0e-14);
  EXPECT_NEAR(nearest_pt[2], 0.0, 1.0e-14);

  // Should project onto face
  pt_to_project = Pt(0.98, 0.37, 20.0);
  nearest_pt = quad.calculateNearestPtOnSurface(pt_to_project);
  EXPECT_NEAR(nearest_pt[0], 0.98, 1.0e-14);
  EXPECT_NEAR(nearest_pt[1], 0.37, 1.0e-14);
  EXPECT_NEAR(nearest_pt[2], 0.0, 1.0e-14);

  // Should project onto face
  pt_to_project = Pt(0.98, 0.37, -20.0);
  nearest_pt = quad.calculateNearestPtOnSurface(pt_to_project);
  EXPECT_NEAR(nearest_pt[0], 0.98, 1.0e-14);
  EXPECT_NEAR(nearest_pt[1], 0.37, 1.0e-14);
  EXPECT_NEAR(nearest_pt[2], 0.0, 1.0e-14);
}
}  // namespace
