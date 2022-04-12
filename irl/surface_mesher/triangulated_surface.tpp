// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_SURFACE_MESH_TRIANGULATED_SURFACE_TPP_
#define IRL_SURFACE_MESH_TRIANGULATED_SURFACE_TPP_

namespace IRL {

inline void TriangulatedSurfaceOutput::addVertex(const Pt& a_vertex) {
  vertices_m.push_back(a_vertex);
}

inline void TriangulatedSurfaceOutput::addTriangle(const UnsignedIndex_t a,
                                                   const UnsignedIndex_t b,
                                                   const UnsignedIndex_t c) {
  triangles_m.push_back(
      ProxyTri<PointStorage>::fromNoExistencePlane(vertices_m, {a, b, c}));
}

inline TriangulatedSurfaceOutput::PointStorage&
TriangulatedSurfaceOutput::getVertexList(void) {
  return vertices_m;
}

inline const TriangulatedSurfaceOutput::PointStorage&
TriangulatedSurfaceOutput::getVertexList(void) const {
  return vertices_m;
}

inline TriangulatedSurfaceOutput::TriangleStorage&
TriangulatedSurfaceOutput::getTriangleList(void) {
  return triangles_m;
}

inline const TriangulatedSurfaceOutput::TriangleStorage&
TriangulatedSurfaceOutput::getTriangleList(void) const {
  return triangles_m;
}

inline TriangulatedSurfaceOutput::PointStorage::size_type
TriangulatedSurfaceOutput::nVertices(void) const {
  return vertices_m.size();
}

inline TriangulatedSurfaceOutput::TriangleStorage::size_type
TriangulatedSurfaceOutput::nTriangles(void) const {
  return triangles_m.size();
}

inline void TriangulatedSurfaceOutput::clearVertices(void) {
  vertices_m.clear();
}
inline void TriangulatedSurfaceOutput::clearTriangles(void) {
  triangles_m.clear();
}
inline void TriangulatedSurfaceOutput::clearAll(void) {
  vertices_m.clear();
  triangles_m.clear();
}

inline void TriangulatedSurfaceOutput::refineSize(
    const double a_max_size, const UnsignedIndex_t a_compute_dim,
    std::function<double(const double a_x, const double a_y)> a_func) {
  const auto original_number_of_tris = this->nTriangles();
  UnsignedIndex_t new_pos = static_cast<UnsignedIndex_t>(this->nVertices());

  std::array<UnsignedIndex_t, 2> dims =
      a_compute_dim == 0   ? std::array<UnsignedIndex_t, 2>{{1, 2}}
      : a_compute_dim == 1 ? std::array<UnsignedIndex_t, 2>{{0, 2}}
                           : std::array<UnsignedIndex_t, 2>{{0, 1}};

  for (std::size_t i = 0; i < this->nTriangles(); ++i) {
    while (triangles_m[i].calculateAbsoluteVolume() > a_max_size) {
      //    if (triangle.calculateAbsoluteVolume() > a_max_size) {
      const auto triangle = triangles_m[i];
      const auto& indices = triangle.getIndexMapping();
      vertices_m.resize(new_pos + 3);

      vertices_m[new_pos] = 0.5 * (triangle[1] - triangle[0]) + triangle[0];
      vertices_m[new_pos + 1] = 0.5 * (triangle[2] - triangle[0]) + triangle[0];
      vertices_m[new_pos + 2] = 0.5 * (triangle[2] - triangle[1]) + triangle[1];

      for (UnsignedIndex_t n = 0; n < 3; ++n) {
        auto& vertex = vertices_m[new_pos + n];
        vertex[a_compute_dim] = a_func(vertex[dims[0]], vertex[dims[1]]);
      }

      triangles_m[i] = ProxyTri<PointStorage>::fromNoExistencePlane(
          vertices_m, {indices[0], new_pos, new_pos + 1});
      triangles_m.push_back(ProxyTri<PointStorage>::fromNoExistencePlane(
          vertices_m, {indices[1], new_pos + 2, new_pos}));
      triangles_m.push_back(ProxyTri<PointStorage>::fromNoExistencePlane(
          vertices_m, {new_pos, new_pos + 2, new_pos + 1}));
      triangles_m.push_back(ProxyTri<PointStorage>::fromNoExistencePlane(
          vertices_m, {indices[2], new_pos + 1, new_pos + 2}));

      new_pos += 3;
    }
  }
}

inline void TriangulatedSurfaceOutput::write(std::string& filename) {
  // binary file
  std::string header_info = filename;
  char head[80];
  std::strncpy(head, header_info.c_str(), sizeof(header_info) - 1);
  char attribute[2] = "0";
  char dummy[4] = "0";
  const unsigned long nTriLong = triangles_m.size();
  std::ofstream myfile;
  myfile.open(filename + ".stl", std::ios::out | std::ios::binary);
  myfile.write(head, sizeof(head));
  myfile.write(reinterpret_cast<const char*>(&nTriLong), 4);

  // write down every triangle
  for (std::size_t i = 0; i < triangles_m.size(); ++i) {
    // normal vector coordinates
    myfile.write(dummy, 4);
    myfile.write(dummy, 4);
    myfile.write(dummy, 4);

    const auto& triangle = triangles_m[i];
    auto points = std::array<float, 9>(
        {static_cast<float>(triangle[0][0]), static_cast<float>(triangle[0][1]),
         static_cast<float>(triangle[0][2]), static_cast<float>(triangle[1][0]),
         static_cast<float>(triangle[1][1]), static_cast<float>(triangle[1][2]),
         static_cast<float>(triangle[2][0]), static_cast<float>(triangle[2][1]),
         static_cast<float>(triangle[2][2])});
    // Write all coordinates
    myfile.write(reinterpret_cast<char*>(points.data()), sizeof(float) * 9);
    myfile.write(attribute, 2);
  }
  myfile.close();
}
}  // namespace IRL

#endif  // IRL_SURFACE_MESH_TRIANGULATED_SURFACE_TPP_
