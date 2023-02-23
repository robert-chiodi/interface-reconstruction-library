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
#include "irl/paraboloid_reconstruction/parametrized_surface.h"

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
    Data<IRL::LocalizedParaboloidLink<double>>*
        a_linked_localized_paraboloids) {
  const BasicMesh& mesh = a_interface.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_linked_localized_paraboloids)(i, j, k) =
            IRL::LocalizedParaboloidLink<double>(&a_cell_localizers(i, j, k),
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
        auto moments = IRL::getNormalizedVolumeMoments<IRL::VolumeMoments>(
            cell, a_interface(i, j, k));
        (*a_liquid_volume_fraction)(i, j, k) =
            moments.volume() / cell.calculateVolume();
        (*a_liquid_centroid)(i, j, k) = moments.centroid();
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
         "OutputDuration", "InterfaceCells");
}

void writeOutDiagnostics(const int a_iteration, const double a_dt,
                         const double a_simulation_time,
                         const Data<double>& a_U, const Data<double>& a_V,
                         const Data<double>& a_W,
                         const Data<double>& a_liquid_volume_fraction,
                         const Data<IRL::Paraboloid>& a_interface,
                         std::chrono::duration<double> a_VOF_duration,
                         std::chrono::duration<double> a_recon_duration,
                         std::chrono::duration<double> a_write_duration) {
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
      "%10d %20.4E %12.3F %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20d"
      "\n",
      a_iteration, a_simulation_time, CFL, liquid_volume_fraction_sum,
      liquid_volume_sum,
      liquid_volume_fraction_sum - initial_liquid_volume_fraction_sum,
      liquid_volume_sum - initial_liquid_volume_sum, a_VOF_duration.count(),
      a_recon_duration.count(), a_write_duration.count(),
      number_of_interface_cells);
}

void writeInterfaceToFile(const Data<double>& a_liquid_volume_fraction,
                          const Data<IRL::Paraboloid>& a_liquid_gas_interface,
                          const double a_time, VTKOutput* a_output) {
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();

  std::vector<IRL::ParametrizedSurfaceOutput> surfaces;
  double radius = 0.25;
  double total_surface = 0.0, avg_mean_curv = 0.0;
  double exact_curv = 2.0 / radius;
  double max_mean_curv_error = 0.0;
  double l2_mean_curv_error = 0.0, l2_counter = 0.0;

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int split_proc = static_cast<int>(std::cbrt(static_cast<double>(size)));
  int imin = 0, nx = mesh.getNx(), jmin = 0, ny = mesh.getNy(), kmin = 0,
      nz = mesh.getNz();

  if (size > 1) {
    if (size == split_proc * split_proc * split_proc) {
      for (int i = 0; i < split_proc; i++) {
        for (int j = 0; j < split_proc; j++) {
          for (int k = 0; k < split_proc; k++) {
            if (i + split_proc * j + split_proc * split_proc * k == rank) {
              imin = i * (mesh.getNx() / split_proc);
              nx = std::min((i + 1) * (mesh.getNx() / split_proc),
                            mesh.getNx()) -
                   imin;
              jmin = j * (mesh.getNy() / split_proc);
              ny = std::min((j + 1) * (mesh.getNy() / split_proc),
                            mesh.getNy()) -
                   jmin;
              kmin = k * (mesh.getNz() / split_proc);
              nz = std::min((k + 1) * (mesh.getNz() / split_proc),
                            mesh.getNz()) -
                   kmin;
            }
          }
        }
      }
    } else {
      imin = rank * (mesh.getNx() / size);
      nx = std::min((rank + 1) * (mesh.getNx() / size), mesh.getNx()) - imin;
      jmin = 0;
      ny = mesh.getNy();
      kmin = 0;
      nz = mesh.getNz();
    }
  }

  for (int i = imin; i < imin + nx; ++i) {
    for (int j = jmin; j < jmin + ny; ++j) {
      for (int k = kmin; k < kmin + nz; ++k) {
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
          if (volume_and_surface.getMoments() > -DBL_MAX) {
            auto surface = volume_and_surface.getSurface();
            double length_scale = std::min(
                std::pow(cell.calculateVolume(), 1.0 / 3.0) / 3.0, 1.0e-2);
            surface.setLengthScale(length_scale);
            if (surface.getSurfaceArea() >
                1.0e-6 * length_scale * length_scale) {
              surfaces.push_back(surface);
            }
            total_surface += surface.getSurfaceArea();
            avg_mean_curv += surface.getMeanCurvatureIntegral();
            max_mean_curv_error = std::max(
                max_mean_curv_error,
                std::abs(surface.getAverageMeanCurvature() - exact_curv));
            // }
            l2_mean_curv_error +=
                (surface.getAverageMeanCurvature() - exact_curv) *
                (surface.getAverageMeanCurvature() - exact_curv);
            l2_counter += 1.0;
          }
        }
      }
    }
  }

  l2_mean_curv_error /= l2_counter;
  l2_mean_curv_error = std::sqrt(l2_mean_curv_error);
  avg_mean_curv /= total_surface;
  total_surface = std::sqrt(total_surface);
  a_output->writeParametrizedInterface(a_time, surfaces);
}
