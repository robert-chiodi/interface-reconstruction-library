// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Austin Han <han.austin@outlook.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This example covers how to create
// a LocalizedSeparatorLink network to
// cut a polyhedron by a plane (as a PlanarSeparator)
//  (localized inside PlanarLocalizers).
// This is done to directly obtain
// volumetric moments for the polyhedron
// internal/external to the separator (below/above the plane)
// for each localizer given.

// This will be shown using a Hexahedron
// that lays across a collection of
// localizers representing RectangularCuboid objects and
// mimicking a uniform Cartesian mesh.
// These PlanarLocalizer planes will be set directly, and in general
// are capable of representing any mesh made up of
// convex polyhedra.

#include <iostream>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/generic_cutting/generic_cutting_definitions.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/moments/separated_volume_moments.h"
#include "irl/moments/volume_moments.h"
#include "irl/planar_reconstruction/planar_localizer.h"
#include "irl/planar_reconstruction/planar_separator.h"

// Helper function that sets up the PlanarLocalizers
// to represent a uniform cubic Cartesian mesh in a
// 3x3x3 neighborhood.
void makeCubicPlanarLocalizer(
    std::array<IRL::PlanarLocalizer, 27>* a_planar_localizer);

// Helper function that links together the individual
// LocalizedSeparatorLink objects to form a mesh
// represented by a bi-direction graph.
void setupLinking(
    std::array<IRL::LocalizedSeparatorLink, 27>* a_localized_separator_link);

// Helper function that returns a linear index for
// the 3x3x3 region with i/j/k in range [0,3]
// and fastest to slowest indices of k/j/i
constexpr IRL::UnsignedIndex_t calculateLinearIndex(
    const IRL::UnsignedIndex_t i, const IRL::UnsignedIndex_t j,
    const IRL::UnsignedIndex_t k) {
  return i * 9 + j * 3 + k;
}

int main(void) {
  std::array<IRL::PlanarLocalizer, 27> planar_localizer;
  std::array<IRL::PlanarSeparator, 27> planar_separator;
  std::array<IRL::LocalizedSeparatorLink, 27> localized_separator_link;

  // Construct the LocalizedSeparatorLink objects,
  // providing the PlanarLocalizer and PlanarSeparator
  // it is comprised of.
  // This is done using pointers, so changes to the
  // objects in the planar_localizer and planar_separator
  // arrays will be reflected in the localized_separator_link.
  for (int i = 0; i < 27; ++i) {
    localized_separator_link[i] =
        IRL::LocalizedSeparatorLink(&planar_localizer[i], &planar_separator[i]);
  }

  // Construct the PlanarLocalizer objects to
  // represent a [-1.5 x 1.5]^3 uniform cubic grid.
  makeCubicPlanarLocalizer(&planar_localizer);

  // Let's now setup the PlanarSeparators involved, which
  // will separate the volume of the hexahedron
  // after it is localized in the region dictated by each PlanarLocalizer
  // in the localized_separator and localized_separator_link.

  // Setup the PlanarSeparator objects for each LocalizedSeparatorLink
  // The added plane will be used to separate the volume localized by
  // each PlanarLocalizer in the associated LocalizedSeparatorLink.
  //
  // For the sake of simplicity,
  // we will specify a flat plane across the middle of
  // the domain with a normal (0.0, 1.0, 0.0) and a distance of 0.0.
  // In general, the PlanarSeparator can have any planar orientation,
  // and no continuity is required between the PlanarSeparator objects in
  // neighboring LocalizedSeparatorLink objects.
  for (int i = 0; i < 27; ++i) {
    planar_separator[i].addPlane(IRL::Plane(IRL::Normal(0.0, 1.0, 0.0), 0.0));
  }

  // Specify the connectivity amongst the different
  // LocalizedSeparatorLink objects. This provides information
  // on which two LocalizedSeparatorLink objects each plane
  // in the contained PlanarLocalizer separates. In a graph
  // sense, this specifies a directed edge from one
  // LocalizedSeparatorLink to another (which can be
  // considered noes in the graph).
  setupLinking(&localized_separator_link);

  // Specify the polyhedron that will be separated, in this case
  // a Hexahedron that intersects all cells in our
  // mimicked mesh.
  IRL::Hexahedron hexahedron;
  hexahedron[0] = IRL::Pt(0.25, -1.00, -1.00);
  hexahedron[1] = IRL::Pt(1.00, 1.00, -1.00);
  hexahedron[2] = IRL::Pt(1.00, 1.00, 1.00);
  hexahedron[3] = IRL::Pt(0.25, -1.00, 1.00);
  hexahedron[4] = IRL::Pt(-0.25, -1.00, -1.00);
  hexahedron[5] = IRL::Pt(-1.00, 1.00, -1.00);
  hexahedron[6] = IRL::Pt(-1.00, 1.00, 1.00);
  hexahedron[7] = IRL::Pt(-0.25, -1.00, 1.00);

  // This routine now performs the intersection of `hexahedron`
  // with the regions created by out LocalizedSeparatorLink objects.
  // Returned from this function will be the VolumeMoments (volume
  // and centroid) of the volume laying in each PlanarLocalizer region
  // separated by the PlanarSeparator paired together in each
  // LocalizedSeparatorLink. The starting location of the cutting,
  // such as localized_separator_link[0] is arbitrary, and just
  // denotes a point in the graph to start from. The volume and
  // centroid can be obtained as
  // phase_moments_LocalizedSeparatorLink[unique_id][phase_id],
  // where unique_id is the ID of the LocalizedSeparatorLink and
  // phase_id is 0 for the volume internal to the PlanarSeparator in
  // that region, or 1 for the volume external to it.
  const auto phase_moments_LocalizedSeparatorLink =
      IRL::getNormalizedVolumeMoments<IRL::TaggedAccumulatedVolumeMoments<
          IRL::SeparatedMoments<IRL::VolumeMoments>>>(
          hexahedron, localized_separator_link[0]);

  std::cout.precision(6);
  std::cout << std::scientific;
  std::cout << '\n';
  std::cout << "Comparison Between Expected and Computed Moments\n";
  std::cout << "================================================\n";
  std::cout << '\n';
  std::cout << "For localizer(1,0,0)\n";
  std::cout << "Internal volume\n";
  std::cout << "  Expected : " << 11.0 / 64.0 << '\n';
  std::cout
      << "  Computed : "
      << phase_moments_LocalizedSeparatorLink[calculateLinearIndex(1, 0, 0)][0]
             .volume()
      << '\n';
  std::cout << "Internal volume centroid\n";
  std::cout << "  Expected : "
            << "( " << 0.0 << ", " << -96.0 / (11.0 * 12.0) << ", " << 0.0
            << " )" << '\n';
  std::cout
      << "  Computed : "
      << phase_moments_LocalizedSeparatorLink[calculateLinearIndex(1, 0, 0)][0]
             .centroid()
      << '\n';
  std::cout << '\n';
  std::cout << "For localizer(1, 1, 1)\n";
  std::cout << "Internal volume\n";
  std::cout << "  Expected : " << 5.0 / 32.0 + 1.0 / 3.0 << '\n';
  std::cout
      << "  Computed : "
      << phase_moments_LocalizedSeparatorLink[calculateLinearIndex(1, 1, 1)][0]
             .volume()
      << '\n';
  std::cout << "Internal volume centroid\n";
  std::cout << "  Expected : "
            << "( " << 0.0 << ", " << -104.0 / 423.0 << ", " << 0.0 << " )"
            << '\n';
  std::cout
      << "  Computed : "
      << phase_moments_LocalizedSeparatorLink[calculateLinearIndex(1, 1, 1)][0]
             .centroid()
      << '\n';

  std::cout << "External volume\n";
  std::cout << "  Expected : " << 0.5 << '\n';
  std::cout
      << "  Computed : "
      << phase_moments_LocalizedSeparatorLink[calculateLinearIndex(1, 1, 1)][1]
             .volume()
      << '\n';
  std::cout << "External volume centroid\n";
  std::cout << "  Expected : "
            << "( " << 0.0 << ", " << 0.25 << ", " << 0.0 << " )" << '\n';
  std::cout
      << "  Computed : "
      << phase_moments_LocalizedSeparatorLink[calculateLinearIndex(1, 1, 1)][1]
             .centroid()
      << '\n';
  std::cout << '\n';
}

// Set the PlanarLocalizers to mimick the 3x3x3 uniform Cartesian mesh.
// Each PlanarLocalizer object will represent a single cell of
// dx = 1.0, with the 3x3x3 region centered at (0.0, 0.0, 0.0)
void makeCubicPlanarLocalizer(
    std::array<IRL::PlanarLocalizer, 27>* a_planar_localizer) {
  for (IRL::UnsignedIndex_t i = 0; i < 3; ++i) {
    for (IRL::UnsignedIndex_t j = 0; j < 3; ++j) {
      for (IRL::UnsignedIndex_t k = 0; k < 3; ++k) {
        // Bounding points for the cube we are creating
        const auto lower_pt = IRL::Pt(static_cast<double>(i) - 1.0 - 0.5,
                                      static_cast<double>(j) - 1.0 - 0.5,
                                      static_cast<double>(k) - 1.0 - 0.5);
        const auto upper_pt = IRL::Pt(static_cast<double>(i) - 1.0 + 0.5,
                                      static_cast<double>(j) - 1.0 + 0.5,
                                      static_cast<double>(k) - 1.0 + 0.5);
        const IRL::UnsignedIndex_t linear_index = calculateLinearIndex(i, j, k);

        // Add the six planes that represent the cube
        (*a_planar_localizer)[linear_index].addPlane(
            IRL::Plane(IRL::Normal(-1.0, 0.0, 0.0), -lower_pt[0]));  // x- face
        (*a_planar_localizer)[linear_index].addPlane(
            IRL::Plane(IRL::Normal(1.0, 0.0, 0.0), upper_pt[0]));  // x+ face

        (*a_planar_localizer)[linear_index].addPlane(
            IRL::Plane(IRL::Normal(0.0, -1.0, 0.0), -lower_pt[1]));  // y- face
        (*a_planar_localizer)[linear_index].addPlane(
            IRL::Plane(IRL::Normal(0.0, 1.0, 0.0), upper_pt[1]));  // y+ face

        (*a_planar_localizer)[linear_index].addPlane(
            IRL::Plane(IRL::Normal(0.0, 0.0, -1.0), -lower_pt[2]));  // z- face
        (*a_planar_localizer)[linear_index].addPlane(
            IRL::Plane(IRL::Normal(0.0, 0.0, 1.0), upper_pt[2]));  // z+ face
      }
    }
  }
}

// Setup the edges between the LocalizedSeparatorLink objects.
// Here, the plane order used in PlanarLocalizer needs to be respected,
// where the PlanarLocalizer objects were setup in the
// makeCubicPlanarLocalizer to be ordered as x-, x+, y-, y+, z-, z+.
// Neighbors therefore neeed to be specified in this way, so the
// 0th neighbor will be the x-1 neighbor. If this neighbor does not
// exist (e.g. for i = 0), a nullptr is supplied instead, which
// is used to indicate a lack of connectivity and termination of the graph.
// A unique id value must also be given to each LocalizedSeparatorLink,
// which is used to tag the moments returned during getNormalizedVolumeMoments
// with the respected region they came from.
void setupLinking(
    std::array<IRL::LocalizedSeparatorLink, 27>* a_localized_separator_link) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        const IRL::UnsignedIndex_t linear_index = calculateLinearIndex(i, j, k);
        (*a_localized_separator_link)[linear_index].setId(linear_index);

        // x-
        IRL::LocalizedSeparatorLink* neighbor =
            i - 1 < 0 ? nullptr
                      : &(*a_localized_separator_link)[calculateLinearIndex(
                            i - 1, j, k)];
        (*a_localized_separator_link)[linear_index].setEdgeConnectivity(
            0, neighbor);

        // x+
        neighbor = i + 1 > 2
                       ? nullptr
                       : &(*a_localized_separator_link)[calculateLinearIndex(
                             i + 1, j, k)];
        (*a_localized_separator_link)[linear_index].setEdgeConnectivity(
            1, neighbor);

        // y-
        neighbor = j - 1 < 0
                       ? nullptr
                       : &(*a_localized_separator_link)[calculateLinearIndex(
                             i, j - 1, k)];
        (*a_localized_separator_link)[linear_index].setEdgeConnectivity(
            2, neighbor);

        // y+
        neighbor = j + 1 > 2
                       ? nullptr
                       : &(*a_localized_separator_link)[calculateLinearIndex(
                             i, j + 1, k)];
        (*a_localized_separator_link)[linear_index].setEdgeConnectivity(
            3, neighbor);

        // z-
        neighbor = k - 1 < 0
                       ? nullptr
                       : &(*a_localized_separator_link)[calculateLinearIndex(
                             i, j, k - 1)];
        (*a_localized_separator_link)[linear_index].setEdgeConnectivity(
            4, neighbor);

        // z+
        neighbor = k + 1 > 2
                       ? nullptr
                       : &(*a_localized_separator_link)[calculateLinearIndex(
                             i, j, k + 1)];
        (*a_localized_separator_link)[linear_index].setEdgeConnectivity(
            5, neighbor);
      }
    }
  }
}
