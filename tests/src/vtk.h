// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_SIMPLE_PARABOLOID_ADVECTOR_VTK_H_
#define EXAMPLES_SIMPLE_PARABOLOID_ADVECTOR_VTK_H_

#include <string>
#include <vector>

#include "irl/paraboloid_reconstruction/parametrized_surface.h"
#include "irl/surface_mesher/triangulated_surface.h"

#include "examples/simple_paraboloid_advector/basic_mesh.h"
#include "examples/simple_paraboloid_advector/data.h"

class VTKOutput {
  struct DataIO {
    DataIO(const std::string& a_name, const Data<double>& a_data)
        : name(a_name), pointer(&a_data) {}

    std::string name;
    const Data<double>* pointer;
  };

 public:
  VTKOutput(const std::string& a_directory, const std::string& a_file_name_base,
            const BasicMesh& a_mesh);

  void addData(const std::string& a_name, const Data<double>& a_data);

  void writeVTKFile(const double a_time);

  void writeVTKInterface(const double a_time,
                         std::vector<IRL::ParametrizedSurfaceOutput>& a_surface,
                         const bool a_print_info = false);

 private:
  std::string directory_m;
  std::string file_name_base_m;
  std::size_t data_files_written_m;
  std::size_t interface_files_written_m;
  const BasicMesh* mesh_m;
  std::vector<DataIO> data_to_write_m;
};

#endif  // EXAMPLES_SIMPLE_PARABOLOID_ADVECTOR_VTK_H_
