// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <stdio.h>
#include <sys/stat.h>
#include <string>

#include "examples/simple_advector/diagnostics.h"

void writeOutDiagnostics(const int a_iteration, const int a_revolution,
                         const int a_step, const int a_steps_per_rev) {
  const double percent_complete = 100.0 * static_cast<double>(a_step) /
                                  static_cast<double>(a_steps_per_rev);
  constexpr int bar_length = 20;
  const int number_of_x =
      static_cast<int>(std::ceil(percent_complete)) / (100 / bar_length);

  printf("%14s %4d %10s", "Revolution: ", a_revolution, "|");
  for (int i = 0; i < number_of_x; ++i) {
    printf("x");
  }
  for (int i = number_of_x; i < bar_length; ++i) {
    printf("-");
  }
  printf("|   %6.2f%%", percent_complete);
  printf("\r");
  fflush(stdout);
}

void newlineDiagnostic(void) { printf("\n"); }

void writeOutMesh(const BasicMesh& a_mesh) {
  std::string output_folder = "viz";
  const int dir_err = mkdir(output_folder.c_str(), 0777);
  FILE* mesh_file;
  std::string file_name = "viz/mesh";
  mesh_file = fopen(file_name.c_str(), "w");
  fprintf(mesh_file, "%10i \n", a_mesh.getNx());
  for (int i = a_mesh.imin(); i <= a_mesh.imax() + 1; ++i) {
    fprintf(mesh_file, "%15.8E \n", a_mesh.x(i));
  }
  fprintf(mesh_file, "%10i \n", a_mesh.getNy());
  for (int j = a_mesh.jmin(); j <= a_mesh.jmax() + 1; ++j) {
    fprintf(mesh_file, "%15.8E \n", a_mesh.y(j));
  }
  fprintf(mesh_file, "%10i \n", a_mesh.getNz());
  for (int k = a_mesh.kmin(); k <= a_mesh.kmax() + 1; ++k) {
    fprintf(mesh_file, "%15.8E \n", a_mesh.z(k));
  }
  fclose(mesh_file);
}

void writeOutVisualization(const int a_viz_output,
                           const Data<double>& a_liquid_volume_fraction) {
  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();
  FILE* viz_file;
  std::string file_name = "viz/vizfile_" + std::to_string(a_viz_output);
  viz_file = fopen(file_name.c_str(), "w");
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        fprintf(viz_file, "%15.8E \n", a_liquid_volume_fraction(i, j, k));
      }
    }
  }
  fclose(viz_file);
}
