// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/paraboloid_advector/solver.h"

#include <stdio.h>
#include <cstdio>
#include <string>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/planar_reconstruction/planar_localizer.h"

#include "examples/paraboloid_advector/basic_mesh.h"
#include "examples/paraboloid_advector/data.h"

// Convert and store the mesh cells into localizers.
void initializeLocalizers(Data<IRL::PlanarLocalizer>* a_localizers) {
  // For each cell in the domain, construct the cell as a RectangularCuboid
  // and get the localizer.
  const BasicMesh& mesh = a_localizers->getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        IRL::Pt lower_corner(mesh.x(i), mesh.y(j), mesh.z(k));
        IRL::Pt upper_corner(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1));
        (*a_localizers)(i, j, k) =
            IRL::RectangularCuboid::fromBoundingPts(lower_corner, upper_corner)
                .getLocalizer();
      }
    }
  }
}

void initializeLocalizedParaboloids(
    const Data<IRL::PlanarLocalizer>& a_cell_localizers,
    const Data<IRL::Paraboloid>& a_interface,
    Data<IRL::LocalizedParaboloidLink>* a_linked_localized_paraboloids) {
  const BasicMesh& mesh = a_interface.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_linked_localized_paraboloids)(i, j, k) =
            IRL::LocalizedParaboloidLink(&a_cell_localizers(i, j, k),
                                         &a_interface(i, j, k));
      }
    }
  }
}

void setPhaseQuantities(const Data<IRL::Paraboloid>& a_interface,
                        Data<double>* a_liquid_volume_fraction,
                        Data<IRL::Pt>* a_liquid_centroid,
                        Data<IRL::Pt>* a_gas_centroid) {
  const BasicMesh& mesh = a_interface.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));
        auto moments = IRL::getNormalizedVolumeMoments<
            IRL::SeparatedMoments<IRL::VolumeMoments>>(cell,
                                                       a_interface(i, j, k));
        (*a_liquid_volume_fraction)(i, j, k) =
            moments[0].volume() / cell.calculateVolume();
        (*a_liquid_centroid)(i, j, k) = moments[0].centroid();
        (*a_gas_centroid)(i, j, k) = moments[1].centroid();
      }
    }
  }
  a_liquid_volume_fraction->updateBorder();
  a_liquid_centroid->updateBorder();
  a_gas_centroid->updateBorder();
}

void writeDiagnosticsHeader(void) {
  printf("%10s %20s %12s %20s %20s %20s %20s %20s %20s %20s\n", "Iteration",
         "Time", "CFL", "liquidVFSum", "liquidVolSum", "ChangeLiquidVFSum",
         "ChangeLiquidVolSum", "AdvectionDuration", "ReconDuration",
         "InterfaceCells");
}

void writeOutDiagnostics(const int a_iteration, const double a_dt,
                         const double a_simulation_time,
                         const Data<double>& a_U, const Data<double>& a_V,
                         const Data<double>& a_W,
                         const Data<double>& a_liquid_volume_fraction,
                         const Data<IRL::Paraboloid>& a_interface,
                         std::chrono::duration<double> a_VOF_duration,
                         std::chrono::duration<double> a_recon_duration) {
  const BasicMesh& mesh = a_U.getMesh();
  static double initial_liquid_volume_fraction_sum;
  static double initial_liquid_volume_sum;
  // Calculate CFL
  double CFL = -DBL_MAX;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        CFL = std::fmax(CFL,
                        std::fmax(a_U(i, j, k) * a_dt / mesh.dx(),
                                  std::fmax(a_V(i, j, k) * a_dt / mesh.dy(),
                                            a_W(i, j, k) * a_dt / mesh.dz())));
      }
    }
  }
  // Calculate sum of volume fraction and sum of liquid volume
  double liquid_volume_fraction_sum = 0.0;
  double liquid_volume_sum = 0.0;
  int number_of_interface_cells = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        liquid_volume_fraction_sum += a_liquid_volume_fraction(i, j, k);
        liquid_volume_sum += a_liquid_volume_fraction(i, j, k) * mesh.dx() *
                             mesh.dy() * mesh.dz();
        if (a_liquid_volume_fraction(i, j, k) >=
                IRL::global_constants::VF_LOW &&
            a_liquid_volume_fraction(i, j, k) <=
                IRL::global_constants::VF_HIGH) {
          ++number_of_interface_cells;
        }
      }
    }
  }
  // Save initial values to compare against.
  if (a_iteration == 0) {
    initial_liquid_volume_fraction_sum = liquid_volume_fraction_sum;
    initial_liquid_volume_sum = liquid_volume_sum;
  }
  printf(
      "%10d %20.4E %12.3F %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20d"
      "\n",
      a_iteration, a_simulation_time, CFL, liquid_volume_fraction_sum,
      liquid_volume_sum,
      liquid_volume_fraction_sum - initial_liquid_volume_fraction_sum,
      liquid_volume_sum - initial_liquid_volume_sum, a_VOF_duration.count(),
      a_recon_duration.count(), number_of_interface_cells);
}

void writeInterfaceToFile(const Data<double>& a_liquid_volume_fraction,
                          const Data<IRL::Paraboloid>& a_liquid_gas_interface,
                          const double a_time, VTKOutput* a_output) {
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();

  std::vector<IRL::TriangulatedSurfaceOutput> surfaces;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        if (a_liquid_volume_fraction(i, j, k) >=
                IRL::global_constants::VF_LOW &&
            a_liquid_volume_fraction(i, j, k) <=
                IRL::global_constants::VF_HIGH) {
          const IRL::Pt lower_cell_pt(mesh.x(i), mesh.y(j), mesh.z(k));
          const IRL::Pt upper_cell_pt(mesh.x(i + 1), mesh.y(j + 1),
                                      mesh.z(k + 1));
          const auto cell = IRL::RectangularCuboid::fromBoundingPts(
              lower_cell_pt, upper_cell_pt);
          auto volume_and_surface = IRL::getVolumeMoments<IRL::AddSurfaceOutput<
              IRL::Volume, IRL::ParametrizedSurfaceOutput>>(
              cell, a_liquid_gas_interface(i, j, k));

          double length_scale =
              std::pow(cell.calculateVolume(), 1.0 / 3.0) / 10.0;
          surfaces.push_back(
              volume_and_surface.getSurface().triangulate(length_scale));
        }
      }
    }
  }
  a_output->writeVTKInterface(a_time, surfaces);
}
