// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/advector/solver.h"

#include <stdio.h>
#include <cstdio>
#include <string>

#include "irl/generic_cutting/generic_cutting.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/polyhedrons/rectangular_cuboid.h"
#include "irl/planar_reconstruction/localized_separator_link.h"
#include "irl/planar_reconstruction/planar_localizer.h"
#include "irl/planar_reconstruction/planar_separator.h"

#include "examples/advector/basic_mesh.h"
#include "examples/advector/data.h"

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

void initializeLocalizedSeparators(
    const Data<IRL::PlanarLocalizer>& a_cell_localizers,
    const Data<IRL::PlanarSeparator>& a_interface,
    Data<IRL::LocalizedSeparatorLink>* a_linked_localized_separators) {
  const BasicMesh& mesh = a_interface.getMesh();
  for (int i = mesh.imino(); i <= mesh.imaxo(); ++i) {
    for (int j = mesh.jmino(); j <= mesh.jmaxo(); ++j) {
      for (int k = mesh.kmino(); k <= mesh.kmaxo(); ++k) {
        (*a_linked_localized_separators)(i, j, k) = IRL::LocalizedSeparatorLink(
            &a_cell_localizers(i, j, k), &a_interface(i, j, k));
      }
    }
  }
}

void setPhaseQuantities(const Data<IRL::PlanarSeparator>& a_interface,
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
  printf("%10s %20s %12s %20s %20s %20s %20s %20s %20s %20s %20s\n",
         "Iteration", "Time", "CFL", "liquidVFSum", "liquidVolSum",
         "ChangeLiquidVFSum", "ChangeLiquidVolSum", "AdvectionDuration",
         "ReconDuration", "InterfaceCells", "TwoPlaneCells");
}

void writeOutMesh(const BasicMesh& a_mesh) {
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

void writeOutDiagnostics(const int a_iteration, const double a_dt,
                         const double a_simulation_time,
                         const Data<double>& a_U, const Data<double>& a_V,
                         const Data<double>& a_W,
                         const Data<double>& a_liquid_volume_fraction,
                         const Data<IRL::PlanarSeparator>& a_interface,
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
  int number_of_two_plane_interface_cells = 0;
  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        liquid_volume_fraction_sum += a_liquid_volume_fraction(i, j, k);
        liquid_volume_sum += a_liquid_volume_fraction(i, j, k) * mesh.dx() *
                             mesh.dy() * mesh.dz();
        if (a_interface(i, j, k).getNumberOfPlanes() > 0 &&
            a_interface(i, j, k)[0].normal().calculateMagnitude() > 0.5) {
          ++number_of_interface_cells;
          if (a_interface(i, j, k).getNumberOfPlanes() > 1) {
            ++number_of_two_plane_interface_cells;
          }
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
      "%10d %20.4E %12.3F %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20d %20d "
      "\n",
      a_iteration, a_simulation_time, CFL, liquid_volume_fraction_sum,
      liquid_volume_sum,
      liquid_volume_fraction_sum - initial_liquid_volume_fraction_sum,
      liquid_volume_sum - initial_liquid_volume_sum, a_VOF_duration.count(),
      a_recon_duration.count(), number_of_interface_cells,
      number_of_two_plane_interface_cells);
}

void writeOutVisualization(const int a_iteration,
                           const int a_visualization_frequency,
                           const double a_simulation_time,
                           const Data<double>& a_liquid_volume_fraction) {
  std::string output_folder = "viz";
  const int dir_err = mkdir(output_folder.c_str(), 0777);

  const BasicMesh& mesh = a_liquid_volume_fraction.getMesh();
  FILE* viz_file;

  // Create zero-filled file name
  int iteration_digits = 6;
  std::string number = std::to_string(a_iteration / a_visualization_frequency);
  assert(iteration_digits > number.length());
  std::string id_suffix =
      std::string(iteration_digits - number.length(), '0') + number;
  std::string file_name = "viz/file_" + id_suffix + ".vtk";
  viz_file = fopen(file_name.c_str(), "w");
  fprintf(viz_file, "# vtk DataFile Version 3.0 \n");
  fprintf(viz_file, "vtk output \n");
  fprintf(viz_file, "ASCII \n");
  fprintf(viz_file, "DATASET RECTILINEAR_GRID \n");
  fprintf(viz_file, "%11s \t %d \t %d \t %d \n", "DIMENSIONS", mesh.getNx() + 1,
          mesh.getNy() + 1, mesh.getNz() + 1);
  fprintf(viz_file, "%14s \t %d \t %6s \n", "X_COORDINATES", mesh.getNx() + 1,
          "float");
  for (int i = mesh.imin(); i <= mesh.imax() + 1; ++i) {
    fprintf(viz_file, "%15.8E \n", mesh.x(i));
  }
  fprintf(viz_file, "%14s \t %d \t %6s \n", "Y_COORDINATES", mesh.getNy() + 1,
          "float");
  for (int j = mesh.jmin(); j <= mesh.jmax() + 1; ++j) {
    fprintf(viz_file, "%15.8E \n", mesh.y(j));
  }
  fprintf(viz_file, "%14s \t %d \t %6s \n", "Z_COORDINATES", mesh.getNz() + 1,
          "float");
  for (int k = mesh.kmin(); k <= mesh.kmax() + 1; ++k) {
    fprintf(viz_file, "%15.8E \n", mesh.z(k));
  }
  fprintf(viz_file, "%11s \t %d \n", "CELL_DATA",
          mesh.getNx() * mesh.getNy() * mesh.getNz());
  fprintf(viz_file, "FIELD FieldData 1 \n");
  fprintf(viz_file, "%11s \t %d \t %d \t %6s \n", "VolumeFraction ", 1,
          mesh.getNx() * mesh.getNy() * mesh.getNz(), "float");
  for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
        fprintf(viz_file, "%15.8E \n", a_liquid_volume_fraction(i, j, k));
      }
    }
  }
  fclose(viz_file);
}

void writeOutInterface(const int a_iteration,
                       const int a_visualization_frequency,
                       const double a_simulation_time,
                       const Data<IRL::PlanarSeparator>& a_interface) {
  std::string output_folder = "viz";
  const int dir_err = mkdir(output_folder.c_str(), 0777);

  const BasicMesh& mesh = a_interface.getMesh();
  FILE* viz_file;

  // Create zero-filled file name
  int iteration_digits = 6;
  std::string number = std::to_string(a_iteration / a_visualization_frequency);
  assert(iteration_digits > number.length());
  std::string id_suffix =
      std::string(iteration_digits - number.length(), '0') + number;
  std::string file_name = "viz/interface_" + id_suffix + ".vtu";
  viz_file = fopen(file_name.c_str(), "w");

  // Build vectors of vertex locations and connectivities
  std::size_t n_vert = 0;
  std::size_t n_faces = 0;
  std::size_t current_face_size = 0;
  std::string vert_loc;
  std::string connectivity;
  std::string offsets;

  auto add_polyhedron = [&](const auto& a_poly) {
    using T = std::decay_t<decltype(a_poly)>;
    std::unordered_map<const typename T::vertex_type*, IRL::UnsignedIndex_t>
        unique_vertices;
    for (IRL::UnsignedIndex_t n = 0; n < a_poly.getNumberOfVertices(); ++n) {
      unique_vertices[a_poly.getVertex(n)] = n + n_vert;
      const auto& vert_pt = a_poly.getVertex(n)->getLocation();
      vert_loc += std::to_string(vert_pt[0]) + " " +
                  std::to_string(vert_pt[1]) + " " +
                  std::to_string(vert_pt[2]) + "\n";
    }
    assert(unique_vertices.size() == a_poly.getNumberOfVertices());

    for (IRL::UnsignedIndex_t n = 0; n < a_poly.getNumberOfFaces(); ++n) {
      const auto& face = a_poly[n];
      auto current_half_edge = face->getStartingHalfEdge();
      do {
        ++current_face_size;
        connectivity +=
            std::to_string(unique_vertices[current_half_edge->getVertex()]) +
            " ";
        current_half_edge = current_half_edge->getNextHalfEdge();
      } while (current_half_edge != face->getStartingHalfEdge());
      offsets += std::to_string(current_face_size) + " ";
      connectivity += "\n";
    }

    n_vert += a_poly.getNumberOfVertices();
    n_faces += a_poly.getNumberOfFaces();
  };

  for (int i = mesh.imin(); i <= mesh.imax(); ++i) {
    for (int j = mesh.jmin(); j <= mesh.jmax(); ++j) {
      for (int k = mesh.kmin(); k <= mesh.kmax(); ++k) {
        const auto& recon = a_interface(i, j, k);
        auto cell = IRL::RectangularCuboid::fromBoundingPts(
            IRL::Pt(mesh.x(i), mesh.y(j), mesh.z(k)),
            IRL::Pt(mesh.x(i + 1), mesh.y(j + 1), mesh.z(k + 1)));

        if (recon.isFlipped()) {
          auto he_poly = cell.generateHalfEdgeVersion();
          auto seg = he_poly.generateSegmentedPolyhedron();
          for (const auto& plane : recon) {
            decltype(seg) clipped;
            auto new_plane = plane.generateFlippedPlane();
            IRL::splitHalfEdgePolytope(&seg, &clipped, &he_poly, new_plane);
            add_polyhedron(clipped);
          }
        } else {
          auto he_poly = cell.generateHalfEdgeVersion();
          auto seg = he_poly.generateSegmentedPolyhedron();
          for (const auto& plane : recon) {
            decltype(seg) clipped;
            IRL::splitHalfEdgePolytope(&seg, &clipped, &he_poly, plane);
          }
          add_polyhedron(seg);
        }
      }
    }
  }

  // Write header
  {
    fprintf(viz_file, "<?xml version=\"1.0\"?>\n");
    fprintf(viz_file,
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n");
    fprintf(viz_file, "<UnstructuredGrid>\n");
    fprintf(viz_file, "<Piece NumberOfPoints=\"%zu\" NumberOfCells=\"%zu\">\n",
            n_vert, n_faces);

    fprintf(viz_file, "<Points>\n");
    fprintf(viz_file,
            "<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n");
    fprintf(viz_file, "%s", vert_loc.c_str());
    fprintf(viz_file, "</DataArray>\n");
    fprintf(viz_file, "</Points>\n");

    fprintf(viz_file, "<Cells>\n");
    fprintf(viz_file,
            "<DataArray type=\"Int32\" Name=\"connectivity\" "
            "format=\"ascii\">\n");
    fprintf(viz_file, "%s", connectivity.c_str());
    fprintf(viz_file, "</DataArray>\n");

    fprintf(viz_file,
            "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    fprintf(viz_file, "%s", offsets.c_str());
    fprintf(viz_file, "\n</DataArray>\n");

    // Cell type - General Polygon type
    fprintf(viz_file,
            "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (std::size_t n = 0; n < n_faces; ++n) {
      fprintf(viz_file, "7 ");  // General polygon type
    }
    fprintf(viz_file,
            "\n</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</"
            "VTKFile>\n");
  }
  fclose(viz_file);
}
