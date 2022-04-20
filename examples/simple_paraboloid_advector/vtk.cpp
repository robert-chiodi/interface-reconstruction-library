
#include "examples/simple_paraboloid_advector/vtk.h"

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
    const double a_time,
    const std::vector<IRL::TriangulatedSurfaceOutput>& a_surface) {
  const auto file_name = directory_m + "/" + file_name_base_m + "_interface_" +
                         std::to_string(interface_files_written_m) + ".vtu";

  FILE* file;
  file = fopen(file_name.c_str(), "w");

  std::size_t number_of_vertices = 0;
  for (std::size_t i = 0; i < a_surface.size(); ++i) {
    const auto& vlist = a_surface[i].getVertexList();
    number_of_vertices += vlist.size();
  }

  std::size_t number_of_triangles = 0;
  for (std::size_t i = 0; i < a_surface.size(); ++i) {
    const auto& tlist = a_surface[i].getTriangleList();
    number_of_triangles += tlist.size();
  }

  fprintf(file, "<VTKFile type=\"UnstructuredGrid\">\n");
  fprintf(file, "<UnstructuredGrid>\n");
  fprintf(file, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
          static_cast<int>(number_of_vertices),
          static_cast<int>(number_of_triangles));
  fprintf(file, "<Points>\n");
  fprintf(file, "<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n");
  std::vector<std::size_t> offset(a_surface.size() + 1, 0);
  for (std::size_t i = 0; i < a_surface.size(); ++i) {
    const auto& vlist = a_surface[i].getVertexList();
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

  for (std::size_t i = 0; i < a_surface.size(); ++i) {
    const auto& tlist = a_surface[i].getTriangleList();
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

  fprintf(file, "</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");

  ++interface_files_written_m;
}
