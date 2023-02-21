// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2023 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "tools/src/surface_io/surface_io.h"

#include <limits.h>
#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

std::vector<std::string> split(const std::string& s, char c);

int main(int argc, char* argv[]) {
  if (argc < 2) {
    simpleErrorHandler("Program expects name of file(s).");
  }

  bool print_info = false;

  IRL::UnsignedIndex_t start_files = 1;
  double length_scale = 0.0;
  if (std::string(argv[1]) == "-l") {
    length_scale = std::atof(argv[2]);
    start_files = 3;
  }
  std::vector<std::pair<int, std::filesystem::path>> files;
  for (IRL::UnsignedIndex_t i = start_files; i < argc; i++) {
    if (std::filesystem::exists(argv[i])) {
      files.push_back({INT_MIN, std::filesystem::path(argv[i])});
    }
  }

  std::filesystem::path input_extension = "irl";
  std::filesystem::path output_extension = "vtu";
  IRL::UnsignedIndex_t number_of_files = files.size();
  for (IRL::UnsignedIndex_t i = 0; i < number_of_files; i++) {
    std::string str = files[i].second.replace_extension("").string();
    std::size_t last_index = str.find_last_not_of("0123456789");
    try {
      files[i].first = std::stoi(str.substr(last_index + 1));
    } catch (...) {
    }
  }
  std::sort(files.begin(), files.end());

  bool check_overwrite = true;
  for (IRL::UnsignedIndex_t i = 0; i < number_of_files; i++) {
    std::string input_file =
        files[i].second.replace_extension(input_extension).string();
    std::string output_file =
        files[i].second.replace_extension(output_extension).string();
    std::string edge_output_file =
        files[i]
            .second
            .replace_filename(std::filesystem::path(
                "edge_" + files[i].second.filename().string()))
            .string();
    if (check_overwrite && std::filesystem::exists(output_file)) {
      std::cout << "  File " << output_file
                << " already exists. Do you want to overwrite? (Y/N) or YALL: ";
      std::string response;
      std::getline(std::cin, response);
      if (response == "YALL") {
        check_overwrite = false;
      } else if (!(response == "Y" || response == "y" || response == "Yes" ||
                   response == "yes")) {
        continue;
      }
    }
    std::cout << "Converting file " << input_file << std::endl;

    /************************ CONVERT *************************/
    std::ifstream stream_input_file(input_file);
    if (stream_input_file) {
      std::string line;
      std::vector<std::string> vec;
      getline(stream_input_file, line);
      for (const std::string& str : split(line, ' ')) {
        vec.push_back(str);
      }
      IRL::UnsignedIndex_t number_of_patches = 0;
      if (vec.size() == 5) {
        try {
          number_of_patches = std::stoi(vec[4]);
        } catch (...) {
          std::cout << "Wrong IRL surface format: could not find number of "
                       "surface patches."
                    << std::endl;
          vec.clear();
          stream_input_file.close();
          exit(-1);
        }
      } else {
        std::cout << "Wrong IRL surface format: could not find number of "
                     "surface patches."
                  << std::endl;
        vec.clear();
        stream_input_file.close();
        exit(-1);
      }
      vec.clear();

      std::vector<IRL::ParametrizedSurfaceOutput> surface_patches;
      surface_patches.resize(number_of_patches);
      for (IRL::UnsignedIndex_t j = 0; j < number_of_patches; j++) {
        /*********** Number of arcs */
        getline(stream_input_file, line);
        for (const std::string& str : split(line, ' ')) {
          vec.push_back(str);
        }
        IRL::UnsignedIndex_t number_of_arcs = 0;
        if (vec.size() == 4) {
          try {
            number_of_arcs = std::stoi(vec[3]);
          } catch (...) {
            std::cout << "Wrong IRL surface format: could not find number of "
                         "arcs of surface patch "
                      << j << std::endl;
            vec.clear();
            stream_input_file.close();
            exit(-1);
          }
        } else {
          std::cout << "Wrong IRL surface format: could not find number of "
                       "arcs of surface patch "
                    << j << std::endl;
          vec.clear();
          stream_input_file.close();
          exit(-1);
        }
        vec.clear();

        /*********** Reference frame */
        getline(stream_input_file, line);
        for (const std::string& str : split(line, ' ')) {
          vec.push_back(str);
        }
        IRL::ReferenceFrame frame;
        if (vec.size() == 17) {
          try {
            frame[0][0] = std::stof(vec[3]);
            frame[0][1] = std::stof(vec[4]);
            frame[0][2] = std::stof(vec[5]);
            frame[1][0] = std::stof(vec[8]);
            frame[1][1] = std::stof(vec[9]);
            frame[1][2] = std::stof(vec[10]);
            frame[2][0] = std::stof(vec[13]);
            frame[2][1] = std::stof(vec[14]);
            frame[2][2] = std::stof(vec[15]);
          } catch (...) {
            std::cout << "Wrong IRL surface format: could not find reference "
                         "frame of surface patch "
                      << j << std::endl;
            vec.clear();
            stream_input_file.close();
            exit(-1);
          }
        } else {
          std::cout << "Wrong IRL surface format: could not find reference "
                       "frame of surface patch "
                    << j << std::endl;
          vec.clear();
          stream_input_file.close();
          exit(-1);
        }
        vec.clear();

        /*********** Datum */
        getline(stream_input_file, line);
        for (const std::string& str : split(line, ' ')) {
          vec.push_back(str);
        }
        IRL::Pt datum;
        if (vec.size() == 6) {
          try {
            datum[0] = std::stof(vec[2]);
            datum[1] = std::stof(vec[3]);
            datum[2] = std::stof(vec[4]);
          } catch (...) {
            std::cout << "Wrong IRL surface format: could not find datum of "
                         "surface patch "
                      << j << std::endl;
            vec.clear();
            stream_input_file.close();
            exit(-1);
          }
        } else {
          std::cout << "Wrong IRL surface format: could not find datum of "
                       "surface patch "
                    << j << std::endl;
          vec.clear();
          stream_input_file.close();
          exit(-1);
        }
        vec.clear();

        /*********** Paraboloid */
        getline(stream_input_file, line);
        for (const std::string& str : split(line, ' ')) {
          vec.push_back(str);
        }
        double coeff_a, coeff_b;
        if (vec.size() == 5) {
          try {
            coeff_a = std::stof(vec[2]);
            coeff_b = std::stof(vec[3]);
          } catch (...) {
            std::cout << "Wrong IRL surface format: could not find aligned "
                         "paraboloid of surface patch "
                      << j << std::endl;
            vec.clear();
            stream_input_file.close();
            exit(-1);
          }
        } else {
          std::cout << "Wrong IRL surface format: could not find aligned "
                       "paraboloid of surface patch "
                    << j << std::endl;
          vec.clear();
          stream_input_file.close();
          exit(-1);
        }
        vec.clear();

        /*********** Arcs */
        surface_patches[j] = IRL::ParametrizedSurfaceOutput(
            IRL::Paraboloid(datum, frame, coeff_a, coeff_b));
        for (IRL::UnsignedIndex_t k = 0; k < number_of_arcs; k++) {
          getline(stream_input_file, line);
          for (const std::string& str : split(line, ' ')) {
            vec.push_back(str);
          }
          int index_start, index_end;
          double weight;
          IRL::Pt start_point, control_point, end_point;
          if (vec.size() == 20) {
            try {
              index_start = std::stoi(vec[2]);
              index_end = std::stoi(vec[3]);
              weight = std::stof(vec[4]);
              start_point = IRL::Pt(std::stof(vec[6]), std::stof(vec[7]),
                                    std::stof(vec[8]));
              control_point = IRL::Pt(std::stof(vec[11]), std::stof(vec[12]),
                                      std::stof(vec[13]));
              end_point = IRL::Pt(std::stof(vec[16]), std::stof(vec[17]),
                                  std::stof(vec[18]));
            } catch (...) {
              std::cout << "Wrong IRL surface format: could not find arc " << k
                        << " of surface patch " << j << std::endl;
              vec.clear();
              stream_input_file.close();
              exit(-1);
            }
          } else {
            std::cout << "Wrong IRL surface format: could not find arc " << k
                      << " of surface patch " << j << std::endl;
            vec.clear();
            stream_input_file.close();
            exit(-1);
          }
          surface_patches[j].addArc(
              IRL::RationalBezierArc(start_point, control_point, end_point,
                                     index_start, index_end, weight));
          vec.clear();
        }
      }
      stream_input_file.close();

      //   std::cout << "Triangulating with length-scale " << length_scale
      //             << std::endl;
      /****************** TRIANGULATE */
      std::vector<IRL::TriangulatedSurfaceOutput> triangulated_surface;
      triangulated_surface.resize(number_of_patches);
      double avg_length = 0.0;
      for (IRL::UnsignedIndex_t j = 0; j < number_of_patches; j++) {
        avg_length += std::fabs(surface_patches[j].getSurfaceArea());
      }
      avg_length /= IRL::safelyEpsilon(static_cast<double>(number_of_patches));
      avg_length = std::sqrt(avg_length) / 3.0;
      for (IRL::UnsignedIndex_t j = 0; j < number_of_patches; j++) {
        if (length_scale <= 0.0) {
          auto curv = std::fabs(surface_patches[j].getAverageMeanCurvature());
          surface_patches[j].setLengthScale(std::min(0.1 / curv, avg_length));
        }
        triangulated_surface[j] =
            surface_patches[j].triangulate(length_scale, 1);
      }

      //   std::cout << "Writing file " << output_file << std::endl;
      /***************** WRITE TO FILE*/
      IRL::UnsignedIndex_t number_of_vertices = 0;
      std::vector<IRL::UnsignedIndex_t> offset(triangulated_surface.size() + 1,
                                               0);
      for (IRL::UnsignedIndex_t j = 0; j < triangulated_surface.size(); ++j) {
        const auto& vlist = triangulated_surface[j].getVertexList();
        number_of_vertices += vlist.size();
        offset[j + 1] = offset[j] + vlist.size();
      }

      IRL::UnsignedIndex_t number_of_triangles = 0;
      for (IRL::UnsignedIndex_t j = 0; j < triangulated_surface.size(); ++j) {
        const auto& tlist = triangulated_surface[j].getTriangleList();
        number_of_triangles += tlist.size();
      }

      IRL::UnsignedIndex_t number_of_bdy_edges = 0;
      std::vector<IRL::UnsignedIndex_t> edge_mapping(number_of_vertices, 0);
      for (int j = 0; j < triangulated_surface.size(); ++j) {
        const auto& elist = triangulated_surface[j].getBoundaryEdgeList();
        const auto off = offset[j];
        for (const auto& edge : elist) {
          edge_mapping[off + edge.first] = 1;
          edge_mapping[off + edge.second] = 1;
        }
        number_of_bdy_edges += elist.size();
      }
      IRL::UnsignedIndex_t number_of_bdy_vertices = 0;
      for (IRL::UnsignedIndex_t j = 0; j < number_of_vertices; ++j) {
        if (edge_mapping[j]) {
          edge_mapping[j] = 1 + number_of_bdy_vertices++;
        }
      }

      FILE* file;
      file = fopen(output_file.c_str(), "w");
      fprintf(file, "<VTKFile type=\"UnstructuredGrid\">\n");
      fprintf(file, "<UnstructuredGrid>\n");
      fprintf(file, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
              static_cast<int>(number_of_vertices),
              static_cast<int>(number_of_triangles));
      fprintf(file, "<Points>\n");
      fprintf(file, "<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n");
      for (IRL::UnsignedIndex_t j = 0; j < triangulated_surface.size(); ++j) {
        const auto& vlist = triangulated_surface[j].getVertexList();
        for (const auto& vertex : vlist) {
          fprintf(file, "%15.8E %15.8E %15.8E ", static_cast<double>(vertex[0]),
                  static_cast<double>(vertex[1]),
                  static_cast<double>(vertex[2]));
        }
      }
      fprintf(file, "</DataArray>\n</Points>\n");
      fprintf(file, "<Cells>\n");
      fprintf(file,
              "<DataArray type=\"Int32\" Name=\"connectivity\" "
              "format=\"ascii\">\n");
      for (IRL::UnsignedIndex_t j = 0; j < triangulated_surface.size(); ++j) {
        const auto& tlist = triangulated_surface[j].getTriangleList();
        const auto off = offset[j];
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
      for (IRL::UnsignedIndex_t j = 0; j < number_of_triangles; ++j) {
        fprintf(file, "%d ", static_cast<int>(3 * (j + 1)));
      }
      fprintf(file, "</DataArray>\n");

      fprintf(file,
              "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
      for (IRL::UnsignedIndex_t j = 0; j < number_of_triangles; ++j) {
        fprintf(file, "5 ");
      }
      fprintf(file, "</DataArray>\n");
      fprintf(file, "</Cells>\n");
      if (print_info) {
        fprintf(file, "<PointData>\n");
        fprintf(file,
                "<DataArray type=\"Float64\" Name=\"Normal\" "
                "NumberOfComponents=\"3\" "
                "format=\"ascii\">\n");
        for (IRL::UnsignedIndex_t j = 0; j < number_of_patches; ++j) {
          const auto& vlist = triangulated_surface[j].getVertexList();
          for (const auto& vertex : vlist) {
            auto normal = surface_patches[j].getNormalNonAligned(vertex);
            fprintf(
                file, "%15.8E %15.8E %15.8E ", static_cast<double>(normal[0]),
                static_cast<double>(normal[1]), static_cast<double>(normal[2]));
          }
        }
        fprintf(file, "\n</DataArray>\n");
        fprintf(file,
                "<DataArray type=\"Float64\" Name=\"MeanCurvature\" "
                "format=\"ascii\">\n");
        for (IRL::UnsignedIndex_t j = 0; j < number_of_patches; ++j) {
          const auto& vlist = triangulated_surface[j].getVertexList();
          for (const auto& vertex : vlist) {
            auto mean_curv =
                surface_patches[j].getMeanCurvatureNonAligned(vertex);
            fprintf(file, "%15.8E ", static_cast<double>(mean_curv));
          }
        }
        fprintf(file, "\n</DataArray>\n");
        fprintf(file,
                "<DataArray type=\"Float64\" Name=\"GaussianCurvature\" "
                "format=\"ascii\">\n");
        for (IRL::UnsignedIndex_t j = 0; j < number_of_patches; ++j) {
          const auto& vlist = triangulated_surface[j].getVertexList();
          for (const auto& vertex : vlist) {
            auto gaussian_curv =
                surface_patches[j].getGaussianCurvatureNonAligned(vertex);
            fprintf(file, "%15.8E ", static_cast<double>(gaussian_curv));
          }
        }
        fprintf(file, "\n</DataArray>\n");
        fprintf(file, "</PointData>\n");
        fprintf(file, "<CellData>\n");
        fprintf(
            file,
            "<DataArray type=\"Int32\" Name=\"SurfaceID\" format=\"ascii\">\n");
        for (IRL::UnsignedIndex_t j = 0; j < number_of_patches; ++j) {
          const auto& tlist = triangulated_surface[j].getTriangleList();
          for (std::size_t k = 0; k < tlist.size(); ++k) {
            fprintf(file, "%d ", static_cast<int>(j));
          }
        }
        fprintf(file, "</DataArray>\n");
        fprintf(file,
                "<DataArray type=\"Float64\" Name=\"SurfaceArea\" "
                "format=\"ascii\">\n");
        for (IRL::UnsignedIndex_t j = 0; j < number_of_patches; ++j) {
          auto surface_area = surface_patches[j].getSurfaceArea();
          const auto& tlist = triangulated_surface[j].getTriangleList();
          for (IRL::UnsignedIndex_t k = 0; k < tlist.size(); ++k) {
            fprintf(file, "%15.8E ", static_cast<double>(surface_area));
          }
        }
        fprintf(file, "</DataArray>\n");
        fprintf(file,
                "<DataArray type=\"Float64\" Name=\"AvgNormal\" "
                "NumberOfComponents=\"3\" format=\"ascii\">\n");
        for (IRL::UnsignedIndex_t j = 0; j < number_of_patches; ++j) {
          auto avg_normal = surface_patches[j].getAverageNormal();
          const auto& datum = surface_patches[j].getParaboloid().getDatum();
          const auto& ref_frame =
              surface_patches[j].getParaboloid().getReferenceFrame();
          const auto base_normal = avg_normal;
          avg_normal = IRL::Normal(0.0, 0.0, 0.0);
          for (std::size_t d = 0; d < 3; ++d) {
            for (std::size_t n = 0; n < 3; ++n) {
              avg_normal[n] += ref_frame[d][n] * base_normal[d];
            }
          }
          const auto& tlist = triangulated_surface[j].getTriangleList();
          for (IRL::UnsignedIndex_t k = 0; k < tlist.size(); ++k) {
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
        for (IRL::UnsignedIndex_t j = 0; j < number_of_patches; ++j) {
          auto avg_mean_curv = surface_patches[j].getAverageMeanCurvature();
          const auto& tlist = triangulated_surface[j].getTriangleList();
          for (IRL::UnsignedIndex_t k = 0; k < tlist.size(); ++k) {
            fprintf(file, "%15.8E ", static_cast<double>(avg_mean_curv));
          }
        }
        fprintf(file, "</DataArray>\n");
        fprintf(file,
                "<DataArray type=\"Float64\" Name=\"AvgGaussianCurvature\" "
                "format=\"ascii\">\n");
        for (IRL::UnsignedIndex_t j = 0; j < number_of_patches; ++j) {
          auto avg_gaussian_curv =
              surface_patches[j].getAverageGaussianCurvature();
          const auto& tlist = triangulated_surface[j].getTriangleList();
          for (IRL::UnsignedIndex_t k = 0; k < tlist.size(); ++k) {
            fprintf(file, "%15.8E ", static_cast<double>(avg_gaussian_curv));
          }
        }
        fprintf(file, "</DataArray>\n");
        fprintf(file, "</CellData>\n");
      }
      fprintf(file, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
      fclose(file);

      //////////////////// WRITING EDGES
      file = fopen(edge_output_file.c_str(), "w");
      fprintf(file, "<VTKFile type=\"UnstructuredGrid\">\n");
      fprintf(file, "<UnstructuredGrid>\n");
      fprintf(file, "<Piece NumberOfPoints=\"%d\" NumberOfCells =\"%d\">\n",
              static_cast<int>(number_of_bdy_vertices),
              static_cast<int>(number_of_bdy_edges));
      fprintf(file, "<Points>\n");
      fprintf(file, "<DataArray type=\"Float64\" NumberOfComponents =\"3\">\n");
      for (std::size_t j = 0; j < triangulated_surface.size(); ++j) {
        const auto& vlist = triangulated_surface[j].getVertexList();
        const auto off = offset[j];
        for (std::size_t k = 0; k < vlist.size(); ++k) {
          if (edge_mapping[off + k]) {
            fprintf(file, "%15.8E %15.8E %15.8E ",
                    static_cast<double>(vlist[k][0]),
                    static_cast<double>(vlist[k][1]),
                    static_cast<double>(vlist[k][2]));
          }
        }
      }
      fprintf(file, "</DataArray>\n</Points>\n");
      fprintf(file, "<Cells>\n");
      fprintf(file,
              "<DataArray type=\"Int32\" Name=\"connectivity\" "
              "format=\"ascii\">\n");
      for (std::size_t j = 0; j < triangulated_surface.size(); ++j) {
        const auto& elist = triangulated_surface[j].getBoundaryEdgeList();
        const auto off = offset[j];
        for (const auto& edge : elist) {
          fprintf(file, "%d %d ",
                  static_cast<int>(edge_mapping[off + edge.first] - 1),
                  static_cast<int>(edge_mapping[off + edge.second] - 1));
        }
      }
      fprintf(file, "</DataArray>\n");
      fprintf(
          file,
          "<DataArray type=\"Int32\" Name=\"offsets\" format =\"ascii\">\n");
      for (std::size_t j = 0; j < number_of_bdy_edges; ++j) {
        fprintf(file, "%d ", static_cast<int>(2 * (j + 1)));
      }
      fprintf(file, "</DataArray>\n");
      fprintf(file,
              "<DataArray type=\"UInt8\" Name=\"types\" format =\"ascii\">\n");
      for (std::size_t j = 0; j < number_of_bdy_edges; ++j) {
        fprintf(file, "3 ");
      }
      fprintf(file, "</DataArray>\n");
      fprintf(file, "</Cells>\n");
      fprintf(file, "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");
      fclose(file);

      offset.clear();
      edge_mapping.clear();
      triangulated_surface.clear();
      surface_patches.clear();
    } else {
      std::cout << "Could not open " << input_file << std::endl;
    }
  }
  std::cout << "-----> Converted " << number_of_files << " files" << std::endl;

  return 0;
}

std::vector<std::string> split(const std::string& s, char c) {
  std::vector<std::string> splitted;

  std::string word;
  for (char ch : s) {
    if ((ch == c) && (!word.empty())) {
      splitted.push_back(word);
      word.clear();
    } else
      word += ch;
  }
  if (!word.empty()) splitted.push_back(word);

  return splitted;
}

static inline void writeLicense(std::ostream& out) {
  out << "// This file is part of the Interface Reconstruction Library (IRL)"
      << std::endl;
  out << "// a library for interface reconstruction and computational geometry "
         "operations"
      << std::endl;
  out << "//" << std::endl;
  out << "// Copyright (C) 2023 Fabien Evrard <fa.evrard@gmail.com>"
      << std::endl;
  out << "//" << std::endl;
  out << "// This Source Code Form is subject to the terms of the Mozilla "
         "Public"
      << std::endl;
  out << "// License, v. 2.0. If a copy of the MPL was not distributed with "
         "this"
      << std::endl;
  out << "// file, You can obtain one at https://mozilla.org/MPL/2.0/.\n"
      << std::endl;
}