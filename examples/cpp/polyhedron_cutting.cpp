// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi  <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// In this example, the ability to obtain volumetric moments
// for polyhedra is demonstrated using IRL's getNormalizedVolumeMoments
// function.

#include <iostream>

#include "src/generic_cutting/generic_cutting.h"
#include "src/generic_cutting/generic_cutting_definitions.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/general/pt.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/moments/separated_volume_moments.h"
#include "src/moments/volume_moments.h"
#include "src/planar_reconstruction/planar_separator.h"

int main(void) {
  // Set up a polyhedron using IRL's GeneralPolyhedron class.
  // Note: Many specific polyhedron types already exist in IRL
  // with some methods being specialized to take advantage of the
  // assumed connectivity knowledge when they are used (such
  // as more efficient calculation of volume for a cube). The general
  // polyhedron type is the most flexible, however, and works
  // with getVolumeMoments and the reconstruction methods.

  // Specify a boundary representation for connectivity of the vertices.
  std::array<std::array<IRL::UnsignedIndex_t, 4>, 6> polyhedron_brep;
  polyhedron_brep[0] = std::array<IRL::UnsignedIndex_t, 4>{{0, 1, 2, 3}};
  polyhedron_brep[1] = std::array<IRL::UnsignedIndex_t, 4>{{7, 6, 5, 4}};
  polyhedron_brep[2] = std::array<IRL::UnsignedIndex_t, 4>{{0, 3, 7, 4}};
  polyhedron_brep[3] = std::array<IRL::UnsignedIndex_t, 4>{{1, 5, 6, 2}};
  polyhedron_brep[4] = std::array<IRL::UnsignedIndex_t, 4>{{2, 6, 7, 3}};
  polyhedron_brep[5] = std::array<IRL::UnsignedIndex_t, 4>{{0, 4, 5, 1}};
  const auto polyhedron_connectivity =
      IRL::PolyhedronConnectivity(polyhedron_brep);

  // Create a GeneralPolyhedron specifying the location of each
  // vertex and the connectivity. The same connectivity could be reused for
  // multiple objects.
  // Note: In this specific case (of a unit-cell), the polyhedron could have
  // been easily generated using
  // IRL::RectangularCuboid::fromBoundingPts(IRL::Pt(-0.5,-0.5,-0,5),
  //                                         IRL::Pt(0.5,0.5,0.5))
  auto vertex_locations = std::array<IRL::Pt, 8>{
      {IRL::Pt(0.5, -0.5, -0.5), IRL::Pt(0.5, 0.5, -0.5),
       IRL::Pt(0.5, 0.5, 0.5), IRL::Pt(0.5, -0.5, 0.5),
       IRL::Pt(-0.5, -0.5, -0.5), IRL::Pt(-0.5, 0.5, -0.5),
       IRL::Pt(-0.5, 0.5, 0.5), IRL::Pt(-0.5, -0.5, 0.5)}};
  auto cube =
      IRL::GeneralPolyhedron(vertex_locations, &polyhedron_connectivity);

  // Define interface reconstruction representing a x-z sheet
  // centered -0.25 from cell center.
  const auto planar_separator = IRL::PlanarSeparator::fromTwoPlanes(
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
  std::cout << std::scientific;
  std::cout << '\n';
  std::cout << "Comparison between expected and computed results\n";
  std::cout << "================================================\n";
  std::cout << "Volume between planes\n";
  std::cout << "  Expected         : " << 0.1 << '\n';
  std::cout << "  Computed         : " << phase_moments[0].volume() << '\n';
  std::cout << "Centroid for volume between planes\n";
  std::cout << "  Expected         : "
            << "( " << 0.0 << ", " << -0.25 << ", " << 0.0 << " )" << '\n';
  std::cout << "  Computed         : " << phase_moments[0].centroid() << '\n';
  std::cout << "Volume outside of planes\n";
  std::cout << "  Expected         : " << 0.9 << '\n';
  std::cout << "  Computed         : " << phase_moments[1].volume() << '\n';
  std::cout << "Centroid for volume outside of planes\n";
  std::cout << "  Expected         : "
            << "( " << 0.0 << ", " << -(0.1 * (-0.25)) / 0.9 << ", " << 0.0
            << " )" << '\n';
  std::cout << "  Computed         : " << phase_moments[1].centroid() << '\n';
}
