
#include "examples/paraboloid_advector/vtk.h"

#include <stdio.h>
#include <sys/stat.h>

VTKOutput::VTKOutput(const std::string& a_directory,
                     const std::string& a_file_name_base,
                     const BasicMesh& a_mesh)
    : directory_m(a_directory),
      file_name_base_m(a_file_name_base),
      data_files_written_m(0),
      interface_files_written_m(0),
      mesh_m(&a_mesh),
      data_to_write_m() {
  const int dir_err = mkdir(directory_m.c_str(), 0777);
}

void VTKOutput::addData(const std::string& a_name, const Data<double>& a_data) {
  data_to_write_m.push_back(DataIO(a_name, a_data));
}

void VTKOutput::writeVTKFile(const double a_time) {
  const auto file_name = directory_m + "/" + file_name_base_m + "_" +
                         std::to_string(data_files_written_m) + ".vtr";

  FILE* file;
  file = fopen(file_name.c_str(), "w");

  fprintf(file, "<VTKFile type=\"RectilinearGrid\">\n");
  fprintf(file, "<RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">\n",
          mesh_m->imin(), mesh_m->imax() + 1, mesh_m->jmin(),
          mesh_m->jmax() + 1, mesh_m->kmin(), mesh_m->kmax() + 1);
  fprintf(file, "<Piece Extent=\"%d %d %d %d %d %d\">\n", mesh_m->imin(),
          mesh_m->imax() + 1, mesh_m->jmin(), mesh_m->jmax() + 1,
          mesh_m->kmin(), mesh_m->kmax() + 1);

  fprintf(file, "<Coordinates>\n");
  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = mesh_m->imin(); i <= mesh_m->imax() + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(mesh_m->x(i)));
  }
  fprintf(file, "\n</DataArray>\n");
  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = mesh_m->jmin(); i <= mesh_m->jmax() + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(mesh_m->y(i)));
  }
  fprintf(file, "\n</DataArray>\n");

  fprintf(file,
          "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
          "format=\"ascii\">\n");
  for (int i = mesh_m->kmin(); i <= mesh_m->kmax() + 1; ++i) {
    fprintf(file, "%15.8E ", static_cast<double>(mesh_m->z(i)));
  }
  fprintf(file, "\n</DataArray>\n");

  fprintf(file, "</Coordinates>\n");

  fprintf(file, "<PointData>\n</PointData>\n");

  fprintf(file, "<CellData Scalars=\"");
  for (auto& data : data_to_write_m) {
    fprintf(file, "%s ", data.name.c_str());
  }
  fprintf(file, "\" >\n");
  for (auto& data : data_to_write_m) {
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\" "
            "format=\"ascii\">\n",
            data.name.c_str());
    for (int k = mesh_m->kmin(); k <= mesh_m->kmax(); ++k) {
      for (int j = mesh_m->jmin(); j <= mesh_m->jmax(); ++j) {
        for (int i = mesh_m->imin(); i <= mesh_m->imax(); ++i) {
          fprintf(file, "%15.8E ",
                  static_cast<double>((*data.pointer)(i, j, k)));
        }
      }
    }
    fprintf(file, "\n</DataArray>\n");
  }
  fprintf(file, "</CellData>\n");
  fprintf(file, "</Piece>\n</RectilinearGrid>\n</VTKFile>\n");
  fclose(file);
  ++data_files_written_m;
}

void VTKOutput::writeVTKInterface(
    const double a_time, std::vector<IRL::ParametrizedSurfaceOutput>& a_surface,
    const bool a_print_info) {
  const auto surface_file_name =
      directory_m + "/" + file_name_base_m + "_interface_" +
      std::to_string(interface_files_written_m) + ".vtu";
  const auto bdy_edge_file_name =
      directory_m + "/" + file_name_base_m + "_interface_boundary_" +
      std::to_string(interface_files_written_m) + ".vtu";

  FILE* file;

  //////////////////// WRITING TRIANGLES
  file = fopen(surface_file_name.c_str(), "w");

  std::vector<IRL::TriangulatedSurfaceOutput> triangulated_surface;
  triangulated_surface.resize(a_surface.size());
  for (std::size_t i = 0; i < a_surface.size(); ++i) {
    triangulated_surface[i] = a_surface[i].triangulate();
  }

  std::size_t number_of_vertices = 0;
  for (std::size_t i = 0; i < triangulated_surface.size(); ++i) {
    const auto& vlist = triangulated_surface[i].getVertexList();
    number_of_vertices += vlist.size();
  }

  std::size_t number_of_triangles = 0;
  for (std::size_t i = 0; i < triangulated_surface.size(); ++i) {
    const auto& tlist = triangulated_surface[i].getTriangleList();
    number_of_triangles += tlist.size();
  }

  fprintf(file, "<VTKFile type=\"UnstructuredGrid\">\n");
  fprintf(file, "<UnstructuredGrid>\n");
  fprintf(file, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
          static_cast<int>(number_of_vertices),
          static_cast<int>(number_of_triangles));
  fprintf(file, "<Points>\n");
  fprintf(file, "<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n");
  std::vector<std::size_t> offset(triangulated_surface.size() + 1, 0);
  for (std::size_t i = 0; i < triangulated_surface.size(); ++i) {
    const auto& vlist = triangulated_surface[i].getVertexList();
    for (const auto& vertex : vlist) {
      fprintf(file, "%15.8E %15.8E %15.8E ", static_cast<double>(vertex[0]),
              static_cast<double>(vertex[1]), static_cast<double>(vertex[2]));
    }
    offset[i + 1] = offset[i] + vlist.size();
  }
  fprintf(file, "</DataArray>\n</Points>\n");

  fprintf(file, "<Cells>\n");
  fprintf(
      file,
      "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");

  for (std::size_t i = 0; i < triangulated_surface.size(); ++i) {
    const auto& tlist = triangulated_surface[i].getTriangleList();
    const auto off = offset[i];
    for (const auto& triangle : tlist) {
      const auto& index_mapping = triangle.getIndexMapping();
      fprintf(file, "%d %d %d ", static_cast<int>(off + index_mapping[0]),
              static_cast<int>(off + index_mapping[1]),
              static_cast<int>(off + index_mapping[2]));
    }
  }
  fprintf(file, "</DataArray>\n");

  fprintf(file,
          "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  for (std::size_t i = 0; i < number_of_triangles; ++i) {
    fprintf(file, "%d ", static_cast<int>(3 * (i + 1)));
  }
  fprintf(file, "</DataArray>\n");

  fprintf(file, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (std::size_t i = 0; i < number_of_triangles; ++i) {
    fprintf(file, "5 ");
  }
  fprintf(file, "</DataArray>\n");

  fprintf(file, "</Cells>\n");

  if (a_print_info) {
    fprintf(file, "<PointData>\n");
    fprintf(
        file,
        "<DataArray type=\"Float64\" Name=\"Normal\" NumberOfComponents=\"3\" "
        "format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      const auto& vlist = triangulated_surface[i].getVertexList();
      for (const auto& vertex : vlist) {
        auto normal = a_surface[i].getNormalNonAligned(vertex);
        fprintf(file, "%15.8E %15.8E %15.8E ", static_cast<double>(normal[0]),
                static_cast<double>(normal[1]), static_cast<double>(normal[2]));
      }
    }
    fprintf(file, "\n</DataArray>\n");
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"MeanCurvature\" "
            "format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      const auto& vlist = triangulated_surface[i].getVertexList();
      for (const auto& vertex : vlist) {
        auto mean_curv = a_surface[i].getMeanCurvatureNonAligned(vertex);
        fprintf(file, "%15.8E ", static_cast<double>(mean_curv));
      }
    }
    fprintf(file, "\n</DataArray>\n");
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"GaussianCurvature\" "
            "format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      const auto& vlist = triangulated_surface[i].getVertexList();
      for (const auto& vertex : vlist) {
        auto gaussian_curv =
            a_surface[i].getGaussianCurvatureNonAligned(vertex);
        fprintf(file, "%15.8E ", static_cast<double>(gaussian_curv));
      }
    }
    fprintf(file, "\n</DataArray>\n");
    fprintf(file, "</PointData>\n");

    fprintf(file, "<CellData>\n");
    fprintf(file,
            "<DataArray type=\"Int32\" Name=\"SurfaceID\" format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      const auto& tlist = triangulated_surface[i].getTriangleList();
      for (std::size_t j = 0; j < tlist.size(); ++j) {
        fprintf(file, "%d ", static_cast<int>(i));
      }
    }
    fprintf(file, "</DataArray>\n");
    fprintf(
        file,
        "<DataArray type=\"Float64\" Name=\"SurfaceArea\" format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      auto surface_area = a_surface[i].getSurfaceArea();
      const auto& tlist = triangulated_surface[i].getTriangleList();
      for (std::size_t j = 0; j < tlist.size(); ++j) {
        fprintf(file, "%15.8E ", static_cast<double>(surface_area));
      }
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"AvgNormal\" "
            "NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      auto avg_normal = a_surface[i].getAverageNormal();
      const auto& datum = a_surface[i].getParaboloid().getDatum();
      const auto& ref_frame = a_surface[i].getParaboloid().getReferenceFrame();
      const auto base_normal = avg_normal;
      avg_normal = IRL::Normal(0.0, 0.0, 0.0);
      for (std::size_t d = 0; d < 3; ++d) {
        for (std::size_t n = 0; n < 3; ++n) {
          avg_normal[n] += ref_frame[d][n] * base_normal[d];
        }
      }
      const auto& tlist = triangulated_surface[i].getTriangleList();
      for (std::size_t j = 0; j < tlist.size(); ++j) {
        fprintf(file, "%15.8E %15.8E %15.8E ",
                static_cast<double>(avg_normal[0]),
                static_cast<double>(avg_normal[1]),
                static_cast<double>(avg_normal[2]));
      }
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"AvgMeanCurvature\" "
            "format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      auto avg_mean_curv = a_surface[i].getAverageMeanCurvature();
      const auto& tlist = triangulated_surface[i].getTriangleList();
      for (std::size_t j = 0; j < tlist.size(); ++j) {
        fprintf(file, "%15.8E ", static_cast<double>(avg_mean_curv));
      }
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"AvgGaussianCurvature\" "
            "format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      auto avg_gaussian_curv = a_surface[i].getAverageGaussianCurvature();
      const auto& tlist = triangulated_surface[i].getTriangleList();
      for (std::size_t j = 0; j < tlist.size(); ++j) {
        fprintf(file, "%15.8E ", static_cast<double>(avg_gaussian_curv));
      }
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "</CellData>\n");
  }
  fprintf(file, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
  fclose(file);

  //////////////////// WRITING EDGES
  file = fopen(bdy_edge_file_name.c_str(), "w");
  std::size_t number_of_bdy_edges = 0;
  std::vector<std::size_t> edge_mapping(number_of_vertices, 0);
  for (std::size_t i = 0; i < triangulated_surface.size(); ++i) {
    const auto& elist = triangulated_surface[i].getBoundaryEdgeList();
    const auto off = offset[i];
    for (const auto& edge : elist) {
      if (!edge_mapping[off + edge.first]) {
        edge_mapping[off + edge.first] = 1;
      }
      if (!edge_mapping[off + edge.second]) {
        edge_mapping[off + edge.second] = 1;
      }
    }
    number_of_bdy_edges += elist.size();
  }
  std::size_t number_of_bdy_vertices = 0;
  for (std::size_t i = 0; i < number_of_vertices; ++i) {
    if (edge_mapping[i]) {
      edge_mapping[i] = 1 + number_of_bdy_vertices++;
    }
  }

  fprintf(file, "<VTKFile type=\"UnstructuredGrid\">\n");
  fprintf(file, "<UnstructuredGrid>\n");
  fprintf(file, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
          static_cast<int>(number_of_bdy_vertices),
          static_cast<int>(number_of_bdy_edges));
  fprintf(file, "<Points>\n");
  fprintf(file, "<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n");
  for (std::size_t i = 0; i < triangulated_surface.size(); ++i) {
    const auto& vlist = triangulated_surface[i].getVertexList();
    const auto off = offset[i];
    for (std::size_t j = 0; j < vlist.size(); ++j) {
      if (edge_mapping[off + j]) {
        fprintf(file, "%15.8E %15.8E %15.8E ", static_cast<double>(vlist[j][0]),
                static_cast<double>(vlist[j][1]),
                static_cast<double>(vlist[j][2]));
      }
    }
  }
  fprintf(file, "</DataArray>\n</Points>\n");

  fprintf(file, "<Cells>\n");
  fprintf(
      file,
      "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");

  for (std::size_t i = 0; i < triangulated_surface.size(); ++i) {
    const auto& elist = triangulated_surface[i].getBoundaryEdgeList();
    const auto off = offset[i];
    for (const auto& edge : elist) {
      fprintf(file, "%d %d ",
              static_cast<int>(edge_mapping[off + edge.first] - 1),
              static_cast<int>(edge_mapping[off + edge.second] - 1));
    }
  }
  fprintf(file, "</DataArray>\n");

  fprintf(file,
          "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  for (std::size_t i = 0; i < number_of_bdy_edges; ++i) {
    fprintf(file, "%d ", static_cast<int>(2 * (i + 1)));
  }
  fprintf(file, "</DataArray>\n");

  fprintf(file, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (std::size_t i = 0; i < number_of_bdy_edges; ++i) {
    fprintf(file, "3 ");
  }
  fprintf(file, "</DataArray>\n");

  fprintf(file, "</Cells>\n");

  if (a_print_info) {
    fprintf(file, "<PointData>\n");
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"Normal\" "
            "NumberOfComponents=\"3\" "
            "format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      const auto& vlist = triangulated_surface[i].getVertexList();
      const auto off = offset[i];
      for (std::size_t j = 0; j < vlist.size(); ++j) {
        if (edge_mapping[off + j]) {
          auto normal = a_surface[i].getNormalNonAligned(vlist[j]);
          fprintf(file, "%15.8E %15.8E %15.8E ", static_cast<double>(normal[0]),
                  static_cast<double>(normal[1]),
                  static_cast<double>(normal[2]));
        }
      }
    }
    fprintf(file, "\n</DataArray>\n");
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"MeanCurvature\" "
            "format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      const auto& vlist = triangulated_surface[i].getVertexList();
      const auto off = offset[i];
      for (std::size_t j = 0; j < vlist.size(); ++j) {
        if (edge_mapping[off + j]) {
          auto mean_curv = a_surface[i].getMeanCurvatureNonAligned(vlist[j]);
          fprintf(file, "%15.8E ", static_cast<double>(mean_curv));
        }
      }
    }
    fprintf(file, "\n</DataArray>\n");
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"GaussianCurvature\" "
            "format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      const auto& vlist = triangulated_surface[i].getVertexList();
      const auto off = offset[i];
      for (std::size_t j = 0; j < vlist.size(); ++j) {
        if (edge_mapping[off + j]) {
          auto gaussian_curv =
              a_surface[i].getGaussianCurvatureNonAligned(vlist[j]);
          fprintf(file, "%15.8E ", static_cast<double>(gaussian_curv));
        }
      }
    }
    fprintf(file, "\n</DataArray>\n");
    fprintf(file, "</PointData>\n");

    fprintf(file, "<CellData>\n");
    fprintf(file,
            "<DataArray type=\"Int32\" Name=\"SurfaceID\" format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      const auto& elist = triangulated_surface[i].getBoundaryEdgeList();
      for (std::size_t j = 0; j < elist.size(); ++j) {
        fprintf(file, "%d ", static_cast<int>(i));
      }
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"SurfaceArea\" "
            "format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      auto surface_area = a_surface[i].getSurfaceArea();
      const auto& elist = triangulated_surface[i].getBoundaryEdgeList();
      for (std::size_t j = 0; j < elist.size(); ++j) {
        fprintf(file, "%15.8E ", static_cast<double>(surface_area));
      }
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"AvgNormal\" "
            "NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      auto avg_normal = a_surface[i].getAverageNormal();
      const auto& datum = a_surface[i].getParaboloid().getDatum();
      const auto& ref_frame = a_surface[i].getParaboloid().getReferenceFrame();
      const auto base_normal = avg_normal;
      avg_normal = IRL::Normal(0.0, 0.0, 0.0);
      for (std::size_t d = 0; d < 3; ++d) {
        for (std::size_t n = 0; n < 3; ++n) {
          avg_normal[n] += ref_frame[d][n] * base_normal[d];
        }
      }
      const auto& elist = triangulated_surface[i].getBoundaryEdgeList();
      for (std::size_t j = 0; j < elist.size(); ++j) {
        fprintf(file, "%15.8E %15.8E %15.8E ",
                static_cast<double>(avg_normal[0]),
                static_cast<double>(avg_normal[1]),
                static_cast<double>(avg_normal[2]));
      }
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"AvgMeanCurvature\" "
            "format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      auto avg_mean_curv = a_surface[i].getAverageMeanCurvature();
      const auto& elist = triangulated_surface[i].getBoundaryEdgeList();
      for (std::size_t j = 0; j < elist.size(); ++j) {
        fprintf(file, "%15.8E ", static_cast<double>(avg_mean_curv));
      }
    }
    fprintf(file, "</DataArray>\n");
    fprintf(file,
            "<DataArray type=\"Float64\" Name=\"AvgGaussianCurvature\" "
            "format=\"ascii\">\n");
    for (std::size_t i = 0; i < a_surface.size(); ++i) {
      auto avg_gaussian_curv = a_surface[i].getAverageGaussianCurvature();
      const auto& elist = triangulated_surface[i].getBoundaryEdgeList();
      for (std::size_t j = 0; j < elist.size(); ++j) {
        fprintf(file, "%15.8E ", static_cast<double>(avg_gaussian_curv));
      }
    }
    fprintf(file, "</DataArray>\n");

    fprintf(file, "</CellData>\n");
  }
  fprintf(file, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
  fclose(file);

  ++interface_files_written_m;
}
