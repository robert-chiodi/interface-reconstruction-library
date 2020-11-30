// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Austin Han <han.austin@outlook.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// In this example, the moment-of-fluid
// reconstruction (MoF) method implemented inside
// IRL will be demonstrated. First, a PlanarSeparator
// will be directly specified and used to obtain
// SeperatedVolumeMoments for a cell. This will
// then be used to perform a MoF reconstruction.
// The plane from the MoF reconstruction
// is compared to the one initially given.
#include <cmath>
#include <iostream>

#include "irl/generic_cutting/generic_cutting_definitions.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/interface_reconstruction_methods/reconstruction_interface.h"
#include "irl/moments/separated_volume_moments.h"
#include "irl/planar_reconstruction/planar_separator.h"

int main(void) {
  // Perform MoF on a unit cubic cell [-0.5, 0.5]^3
  const IRL::RectangularCuboid cuboid = IRL::unit_cell;

  // Define interface reconstruction representing a slanted
  // line across the cell
  const auto correct_planar_separator = IRL::PlanarSeparator::fromOnePlane(
      IRL::Plane(IRL::Normal::normalized(std::sqrt(3.0), 1.0, 0.0), -0.15));

  // Use IRL's getNormalizedVolumeMoments to calculate the volume
  // and centroid for this plane applied to cuboid
  const auto phase_moments = IRL::getNormalizedVolumeMoments<
      IRL::SeparatedMoments<IRL::VolumeMoments>>(cuboid,
                                                 correct_planar_separator);

  // Now perform a MoF reconstruction and store the result
  // By default, equal weights will be applied in the optimization
  // when trying to recover the phase centroids.
  const auto found_planar_separator_2D_equal =
      reconstructionWithMOF2D(cuboid, phase_moments);
  const auto found_planar_separator_3D_equal =
      reconstructionWithMOF3D(cuboid, phase_moments);

  // We can also supply our own weights for MoF to use, allowing us
  // to bias matching one centroid over the other. Here, let's just
  // try to match the internal centroid.
  const auto found_planar_separator_2D_weighted =
      reconstructionWithMOF2D(cuboid, phase_moments, 1.0, 0.0);
  const auto found_planar_separator_3D_weighted =
      reconstructionWithMOF3D(cuboid, phase_moments, 1.0, 0.0);

  std::cout.precision(6);
  std::cout << std::scientific;
  std::cout << '\n';
  std::cout << "Comparison between given and computed result\n";
  std::cout << "============================================\n";
  std::cout << "Default (equal) weights: \n";
  std::cout << "  Given plane   : " << correct_planar_separator[0];
  std::cout << "  Computed (2D) : " << found_planar_separator_2D_equal[0];
  std::cout << "  Computed (3D) : " << found_planar_separator_3D_equal[0]
            << '\n';
  std::cout << "Set (unequal) weights: \n";
  std::cout << "  Given plane   : " << correct_planar_separator[0];
  std::cout << "  Computed (2D) : " << found_planar_separator_2D_weighted[0];
  std::cout << "  Computed (3D) : " << found_planar_separator_3D_weighted[0]
            << '\n';
}
