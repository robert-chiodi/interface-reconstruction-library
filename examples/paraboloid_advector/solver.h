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

#include "irl/geometry/general/pt.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/planar_reconstruction/planar_localizer.h"

#include "examples/paraboloid_advector/basic_mesh.h"
#include "examples/paraboloid_advector/data.h"
#include "examples/paraboloid_advector/vtk.h"

struct IRLReconstructionsPack {
  Data<IRL::Paraboloid> liquid_gas_interface;
  Data<IRL::PlanarLocalizer> cell_localizers;
  Data<IRL::LocalizedParaboloidLink> traversing_links;
};

void runSimulation(const int a_number_of_cells, const int a_number_of_rotations,
                   const int a_viz_out_freq);

BasicMesh initializeMesh(const int a_number_of_cells);

Data<IRL::Paraboloid> initializeInterface(const BasicMesh& a_mesh);

Data<IRL::PlanarLocalizer> initializeCellLocalizers(const BasicMesh& a_mesh);

Data<IRL::LocalizedParaboloidLink> setupLinkStructure(
    const Data<IRL::Paraboloid>& a_interface,
    const Data<IRL::PlanarLocalizer>& a_cell_localizers);

Data<double> initializeVolumeFraction(const Data<IRL::Paraboloid>& a_interface);

void advanceState(const double a_dt, Data<double>* a_liquid_volume_fraction,
                  IRLReconstructionsPack* a_reconstructions);

void updateVOF(const double a_dt,
               const Data<IRL::LocalizedParaboloidLink>& a_traversing_link,
               Data<double>* a_liquid_volume_fraction);

IRL::Pt backProjectVertex(const double a_dt, const IRL::Pt& a_pt);

void updateReconstruction(const double a_dt,
                          const Data<double>& a_liquid_volume_fraction,
                          Data<IRL::Paraboloid>* a_liquid_gas_interface);

void correctInterfacePlaneBorders(
    Data<IRL::Paraboloid>* a_liquid_gas_interface);

void writeInterfaceToFile(const Data<double>& a_liquid_volume_fraction,
                          const Data<IRL::Paraboloid>& a_liquid_gas_interface,
                          const double a_time, VTKOutput* a_output);

#endif  // EXAMPLES_PARABOLOID_ADVECTOR_SOLVER_H_
