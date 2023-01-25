// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_PARABOLOID_ADVECTOR_SOLVER_H_
#define EXAMPLES_PARABOLOID_ADVECTOR_SOLVER_H_

#include <mpi.h>
#include <sys/stat.h>
#include <chrono>
#include <iostream>
#include <string>

#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/planar_reconstruction/planar_localizer.h"

#include "examples/paraboloid_advector/basic_mesh.h"
#include "examples/paraboloid_advector/data.h"
#include "examples/paraboloid_advector/reconstruction_types.h"
#include "examples/paraboloid_advector/vof_advection.h"
#include "examples/paraboloid_advector/vtk.h"

/// \brief Handles running and advancing the solution according to provided
/// static functions in structs.
template <class SimulationType>
int runSimulation(const std::string& a_advection_method,
                  const std::string& a_reconstruction_method, const double a_dt,
                  const double a_end_time, const int a_visualization_frequency);

// \brief Convert and store the mesh cells into localizers.
void initializeLocalizers(Data<IRL::PlanarLocalizer>* a_localizers);

/// \brief Initialize the linked localized paraboloids used during advection.
void initializeLocalizedParaboloids(
    const Data<IRL::PlanarLocalizer>& a_cell_localizers,
    const Data<IRL::Paraboloid>& a_interface,
    Data<IRL::LocalizedParaboloidLink<double>>* a_linked_localized_paraboloid);

/// \brief Set phase quantities according to the given
void setPhaseQuantities(const Data<IRL::Paraboloid>& a_interface,
                        Data<double>* a_liquid_volume_fraction,
                        Data<IRL::Pt>* a_liquid_centroid,
                        Data<IRL::Pt>* a_gas_centroid);

/// \brief Write out the header for the diagnostics.
void writeDiagnosticsHeader(void);

/// \brief Write out diagnostics to the screen.
void writeOutDiagnostics(const int a_iteration, const double a_dt,
                         const double a_simulation_time,
                         const Data<double>& a_U, const Data<double>& a_V,
                         const Data<double>& a_W,
                         const Data<double>& a_liquid_volume_fraction,
                         const Data<IRL::Paraboloid>& a_interface,
                         std::chrono::duration<double> a_VOF_duration,
                         std::chrono::duration<double> a_recon_duration,
                         std::chrono::duration<double> a_write_duration);

/// \brief Generates triangulated surface and writes to provided VTK file
void writeInterfaceToFile(const Data<double>& a_liquid_volume_fraction,
                          const Data<IRL::Paraboloid>& a_liquid_gas_interface,
                          const double a_time, VTKOutput* a_output);

//******************************************************************* //
//     Template function definitions placed below this.
//******************************************************************* //
template <class SimulationType>
int runSimulation(const std::string& a_advection_method,
                  const std::string& a_reconstruction_method, const double a_dt,
                  const double a_end_time,
                  const int a_visualization_frequency) {
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Set mesh
  BasicMesh cc_mesh = SimulationType::setMesh();

  // Allocate local data
  Data<double> velU(&cc_mesh);
  Data<double> velV(&cc_mesh);
  Data<double> velW(&cc_mesh);
  Data<double> liquid_volume_fraction(&cc_mesh);
  Data<IRL::Pt> liquid_centroid(&cc_mesh);
  Data<IRL::Pt> gas_centroid(&cc_mesh);

  // Allocate
  Data<IRL::PlanarLocalizer> cell_localizers(&cc_mesh);
  initializeLocalizers(&cell_localizers);
  Data<IRL::Paraboloid> interface(&cc_mesh);
  Data<IRL::LocalizedParaboloidLink<double>> link_localized_paraboloids(
      &cc_mesh);
  initializeLocalizedParaboloids(cell_localizers, interface,
                                 &link_localized_paraboloids);
  connectMesh(cc_mesh, &link_localized_paraboloids);
  // Set constants in IRL
  IRL::setMinimumVolumeToTrack(10.0 * DBL_EPSILON * cc_mesh.dx() *
                               cc_mesh.dy() * cc_mesh.dz());
  IRL::setVolumeFractionBounds(1.0e-8);
  IRL::setVolumeFractionTolerance(1.0e-13);

  // Initialize data
  SimulationType::initialize(&velU, &velV, &velW, &interface);
  setPhaseQuantities(interface, &liquid_volume_fraction, &liquid_centroid,
                     &gas_centroid);
  const auto starting_liquid_volume_fraction = liquid_volume_fraction;

  VTKOutput vtk_io("viz_out", "viz", cc_mesh);
  vtk_io.addData("VOF", liquid_volume_fraction);
  vtk_io.addData("VelocityX", velU);
  vtk_io.addData("VelocityY", velV);
  vtk_io.addData("VelocityZ", velW);
  double simulation_time = 0.0;
  int iteration = 0;

  if (rank == 0) {
    vtk_io.writeVTKFile(simulation_time);
  }
  getReconstruction(a_reconstruction_method, liquid_volume_fraction,
                    liquid_centroid, gas_centroid, link_localized_paraboloids,
                    0.0, velU, velV, velW, &interface);

  writeInterfaceToFile(liquid_volume_fraction, interface, simulation_time,
                       &vtk_io);
  if (rank == 0) {
    writeDiagnosticsHeader();
  }
  std::string output_folder = "viz";
  const int dir_err = mkdir(output_folder.c_str(), 0777);
  std::chrono::duration<double> advect_VOF_time(0.0);
  std::chrono::duration<double> recon_time(0.0);
  std::chrono::duration<double> write_time(0.0);
  if (rank == 0) {
    writeOutDiagnostics(iteration, a_dt, simulation_time, velU, velV, velW,
                        liquid_volume_fraction, interface, advect_VOF_time,
                        recon_time, write_time);
  }
  while (simulation_time < a_end_time) {
    const double time_step_to_use =
        std::fmin(a_dt, a_end_time - simulation_time);
    SimulationType::setVelocity(simulation_time + 0.5 * time_step_to_use, &velU,
                                &velV, &velW);

    auto start = std::chrono::system_clock::now();
    advectVOF(a_advection_method, a_reconstruction_method, time_step_to_use,
              velU, velV, velW, &link_localized_paraboloids,
              &liquid_volume_fraction, &liquid_centroid, &gas_centroid,
              &interface);
    auto advect_end = std::chrono::system_clock::now();
    advect_VOF_time = advect_end - start;
    getReconstruction(a_reconstruction_method, liquid_volume_fraction,
                      liquid_centroid, gas_centroid, link_localized_paraboloids,
                      time_step_to_use, velU, velV, velW, &interface);
    auto recon_end = std::chrono::system_clock::now();
    recon_time = recon_end - advect_end;

    if (a_visualization_frequency > 0 &&
        iteration % a_visualization_frequency == 0) {
      if (rank == 0) {
        vtk_io.writeVTKFile(simulation_time);
      }
      writeInterfaceToFile(liquid_volume_fraction, interface, simulation_time,
                           &vtk_io);
    }
    auto write_end = std::chrono::system_clock::now();
    write_time = write_end - recon_end;

    simulation_time += time_step_to_use;
    if (rank == 0) {
      writeOutDiagnostics(iteration + 1, time_step_to_use, simulation_time,
                          velU, velV, velW, liquid_volume_fraction, interface,
                          advect_VOF_time, recon_time, write_time);
    }
    ++iteration;
  }

  // L1 Difference between Starting VOF and ending VOF
  double l1_error = 0.0;
  const BasicMesh& mesh = cc_mesh;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        l1_error += std::fabs(liquid_volume_fraction(i, j, k) -
                              starting_liquid_volume_fraction(i, j, k));
      }
    }
  }
  l1_error /= (static_cast<double>(mesh.getNx() * mesh.getNy() * mesh.getNz()));
  std::cout << "L1 Difference between start and end: " << l1_error << std::endl;

  return 0;
}

#endif  // EXAMPLES_PARABOLOID_ADVECTOR_SOLVER_H_
