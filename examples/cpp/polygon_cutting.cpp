// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Austin Han <han.austin@outlook.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This example demonstrates the ability
// to create a polygon with N-vertices and then
// return individual triangles from its surface.
// The integration of surface moments for localized
// regions is also demonstrated.

#include <iostream>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/generic_cutting_definitions.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polygons/polygon.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/moments/separated_volume_moments.h"
#include "irl/moments/volume_moments.h"
#include "irl/planar_reconstruction/planar_localizer.h"
#include "irl/planar_reconstruction/planar_separator.h"

int main(void) {
  // Set the polygon to be a pentagon in the XZ Plane
  // NOTE: There is an assumed connectivity between
  // subsequent vertices and between the last vertex
  // and the first one
  IRL::Polygon polygon;
  polygon.setNumberOfVertices(5);
  polygon[0] = IRL::Pt(0.0, 0.0, 0.0);
  polygon[1] = IRL::Pt(0.0, 0.0, 1.0);
  polygon[2] = IRL::Pt(0.5, 0.0, 1.5);
  polygon[3] = IRL::Pt(1.0, 0.0, 1.0);
  polygon[4] = IRL::Pt(1.0, 0.0, 0.0);

  // Calculate the plane the polygon exists on
  // Following the right-hand rule, this should be
  // n = (0.0, 1.0, 0.0) and d = 0.0
  polygon.calculateAndSetPlaneOfExistence();
  const auto existence_plane = polygon.getPlaneOfExistence();

  // Surface moments for this polygon can be computed
  // directly or the polygon can be used in
  // getVolumeMoments to first perform
  // truncation by a set of planes.
  const auto full_area = polygon.calculateVolume();

  // Localizer that removes all but the top part of the pentagon.
  const auto planar_localizer = IRL::PlanarLocalizer::fromOnePlane(
      IRL::Plane(IRL::Normal(0.0, 0.0, -1.0), -1.0));
  const auto localized_area =
      IRL::getVolumeMoments<IRL::Volume>(polygon, planar_localizer);

  // Individual triangles can also be obtained from polygon as well
  double triangle_area = {0.0};
  for (IRL::UnsignedIndex_t i = 0;
       i < polygon.getNumberOfSimplicesInDecomposition(); i++) {
    const auto triangle = polygon.getSimplexFromDecomposition(i);
    triangle_area += triangle.calculateVolume();
  }

  // These triangles can also be localized using the PlanarLocalizer object
  double triangle_localized_area = {0.0};
  for (IRL::UnsignedIndex_t i = 0;
       i < polygon.getNumberOfSimplicesInDecomposition(); i++) {
    auto triangle = polygon.getSimplexFromDecomposition(i);
    triangle_localized_area +=
        IRL::getVolumeMoments<IRL::Volume>(triangle, planar_localizer);
  }

  std::cout.precision(6);
  std::cout << std::scientific;
  std::cout << '\n';
  std::cout << "           Polygon usage demonstration           \n";
  std::cout << "=================================================\n";
  std::cout << " Expected existence plane   : "
            << "( 0, 1, 0, 0 )\n";
  std::cout << " Calculated existence plane : " << existence_plane << '\n';
  std::cout << " Correct surface area   : " << 1.0 + 2.0 / 8.0 << '\n';
  std::cout << " Calculated surface area : " << full_area << '\n';
  std::cout << " Calculated surface area from triangles : " << triangle_area
            << '\n';
  std::cout << '\n';
  std::cout << " Correct localized surface area                   : "
            << 2.0 / 8.0 << '\n';
  std::cout << " Calculated localized surface area                : "
            << localized_area << '\n';
  std::cout << " Calculated localized surface area from triangles : "
            << triangle_localized_area << '\n';
}
