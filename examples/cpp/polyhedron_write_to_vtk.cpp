// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi  <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// In this example, the ability to write out polytopes in
// the VTK format is demonstrated for a GeneralPolyhedron
// (set to be a unit cube), that is intersected by two planes.
// The portions clipped by each of the two planes, and the
// remaining unclipped portion, are written in VTK format to
// three separate text files. These files can
// then be visualized with software such as
// VisIt or Paraview.

#include <fstream>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/generic_cutting_definitions.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/half_edge_structures/half_edge_polyhedron.h"
#include "irl/geometry/half_edge_structures/segmented_half_edge_polyhedron.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/moments/separated_volume_moments.h"
#include "irl/moments/volume_moments.h"
#include "irl/planar_reconstruction/planar_separator.h"

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

  // The cube GeneralPolyhedron object can now be
  // split by the planes in the planar_separator,
  // with each written to a VTK file for visualization.
  auto half_edge_cube = cube.generateHalfEdgeVersion();
  auto unclipped_polyhedron = half_edge_cube.generateSegmentedPolyhedron();
  std::vector<decltype(unclipped_polyhedron)> clipped_polyhedrons(
      planar_separator.getNumberOfPlanes());
  for (IRL::UnsignedIndex_t n = 0; n < planar_separator.getNumberOfPlanes();
       ++n) {
    splitHalfEdgePolytope(&unclipped_polyhedron, &(clipped_polyhedrons[n]),
                          &half_edge_cube, planar_separator[n]);
  }
  std::ofstream first_clipped("first_clipped.vtu");
  first_clipped << clipped_polyhedrons[0];
  std::ofstream second_clipped("second_clipped.vtu");
  second_clipped << clipped_polyhedrons[1];
  std::ofstream unclipped("unclipped.vtu");
  unclipped << unclipped_polyhedron;
}
