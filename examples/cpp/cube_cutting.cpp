// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Austin Han <han.austin@outlook.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// In this example a RectangularCuboid
// is cut by a two-plane PlanarSeparator that
// divides the RectangularCuboid into volumes internal
// and external to the separator.
// The computed results compared to the correct
// results are then printed to screen for comparison.
// For use of a GeneralPolyhedron, see the example
// polyhedron_cutting.cpp.

#include <iostream>

#include "src/generic_cutting/generic_cutting.h"
#include "src/generic_cutting/generic_cutting_definitions.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/general/pt.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/moments/separated_volume_moments.h"
#include "src/moments/volume_moments.h"
#include "src/planar_reconstruction/planar_localizer.h"
#include "src/planar_reconstruction/planar_separator.h"

int main(void) {
  // Define unit-cubic cell using IRL RectangularCuboid class
  const auto cube = IRL::RectangularCuboid::fromBoundingPts(
      IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, 0.5));

  // Define interface reconstruction representing a x-z sheet
  // centered -0.25 from the cell center
  auto planar_separator = IRL::PlanarSeparator::fromTwoPlanes(
      IRL::Plane(IRL::Normal(0.0, 1.0, 0.0), -0.2),
      IRL::Plane(IRL::Normal(0.0, -1.0, 0.0), 0.3), 1.0);

  // getNormalizedVolumeMoments is the main work horse of IRL.
  // In general, it is capable of intersecting a
  // polygon or polyhedron by a set of planes to return
  // moments of the intersection between the polygon/polyhedron
  // and the half-spaces of the planes.
  auto phase_moments = IRL::getNormalizedVolumeMoments<
      IRL::SeparatedMoments<IRL::VolumeMoments>>(cube, planar_separator);

  // Print out the computed results.
  // The phase moments in a SeparatedMoments<> object
  // are stored as internal to the PlanarSeparator [0]
  // and external to the PlanarSeparator [1]
  std::cout.precision(6);
  std::cout << std::scientific;
  std::cout << '\n';
  std::cout << "Comparison between expected and computed results\n";
  std::cout << "================================================\n";
  std::cout << "Volume between planes\n";
  std::cout << "  Expected : " << 0.1 << '\n';
  std::cout << "  Computed : " << phase_moments[0].volume() << '\n';
  std::cout << "Centroid for volume between planes\n";
  std::cout << "  Expected : "
            << "( 0, -0.25, 0 )" << '\n';
  std::cout << "  Computed : " << phase_moments[0].centroid() << '\n';
  std::cout << "Volume outside of planes\n";
  std::cout << "  Expected : " << 0.9 << '\n';
  std::cout << "  Computed : " << phase_moments[1].volume() << '\n';
  std::cout << "Centroid for volume outside of planes\n";
  std::cout << "  Expected : "
            << "( 0, " << -(0.1 * (-0.25)) / 0.9 << ", 0 )" << '\n';
  std::cout << "  Computed : " << phase_moments[1].centroid() << '\n';
  std::cout << '\n';
}
