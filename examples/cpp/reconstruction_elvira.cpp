// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Austin Han <han.austin@outlook.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// In this example, the ELVIRA
// reconstruction method implemented inside
// IRL will be demonstrated. First, a PlanarSeparator
// will be directly specified and used to obtain
// SeperatedVolumeMoments for a 3x3 stencil. This will
// then be used to perform a ELVIRA reconstruction
// using the liquid volume fraction. The calculated interface
// is then shown to directly recover the planar interface
// that created the original volume fraction field.
// While this is shown for 2D, the exact same process
// can be done for 3D, with 27 members in the
// ELVIRANeighborhood due to the increase in
// stencil size.
#include <cmath>
#include <iostream>

#include "irl/generic_cutting/generic_cutting_definitions.h"
#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/interface_reconstruction_methods/elvira_neighborhood.h"
#include "irl/interface_reconstruction_methods/reconstruction_interface.h"
#include "irl/planar_reconstruction/planar_separator.h"

int main(void) {
  const auto correct_planar_separator = IRL::PlanarSeparator::fromOnePlane(
      IRL::Plane(IRL::Normal::normalized(std::sqrt(3.0), 1.0, 0.0), -0.15));

  IRL::ELVIRANeighborhood neighborhood;
  neighborhood.resize(9);
  std::array<IRL::RectangularCuboid, 9> cells;
  std::array<double, 9> liquid_volume_fraction;

  // Fill in the ELVIRANeighborhood object.
  // The cells and liquid_volume_fraction are held as
  // pointers, so any change to them after setMember
  // will still be reflected in neighborhood.
  for (int j = -1; j < 2; ++j) {
    for (int i = -1; i < 2; ++i) {
      const double x_shift = static_cast<double>(i);
      const double y_shift = static_cast<double>(j);
      const IRL::UnsignedIndex_t linear_index = (j + 1) * 3 + (i + 1);

      cells[linear_index] = IRL::RectangularCuboid::fromBoundingPts(
          IRL::Pt(-0.5 + x_shift, -0.5 + y_shift, -0.5),
          IRL::Pt(0.5 + x_shift, 0.5 + y_shift, 0.5));
      liquid_volume_fraction[linear_index] =
          getVolumeFraction(cells[linear_index], correct_planar_separator);
      neighborhood.setMember(&cells[linear_index],
                             &liquid_volume_fraction[linear_index], i, j);
    }
  }

  // Now perform actual ELVIRA and obtain interface PlanarSeparator
  auto found_planar_separator = reconstructionWithELVIRA2D(neighborhood);

  std::cout.precision(6);
  std::cout << std::scientific;
  std::cout << '\n';
  std::cout << "Comparison between given and computed normal\n";
  std::cout << "============================================\n";
  std::cout << "Normal vectors\n";
  std::cout << "  Given    : " << correct_planar_separator[0].normal() << '\n';
  std::cout << "  Computed : " << found_planar_separator[0].normal() << "\n\n";
  std::cout << "Plane distance\n";
  std::cout << "  Given    : " << -0.15 << '\n';
  std::cout << "  Computed : " << found_planar_separator[0].distance() << '\n';
  std::cout << '\n';
}
