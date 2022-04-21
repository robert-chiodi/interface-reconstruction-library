// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_SURFACE_MESH_TRIANGULATED_SURFACE_H_
#define IRL_SURFACE_MESH_TRIANGULATED_SURFACE_H_

#include <fstream>
#include <functional>
#include <vector>

#include "irl/geometry/polygons/tri.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

class TriangulatedSurfaceOutput {
 public:
  class PointStorage : public std::vector<Pt> {
   public:
    using pt_type = Pt;
  };
  using EdgeStorage = std::vector<std::pair<UnsignedIndex_t, UnsignedIndex_t>>;
  using TriangleStorage = std::vector<ProxyTri<PointStorage>>;

  /// \brief Default constructor.
  TriangulatedSurfaceOutput(void) = default;
  ~TriangulatedSurfaceOutput(void) = default;

  void addVertex(const Pt& a_vertex);
  void addBoundaryEdge(const UnsignedIndex_t a, const UnsignedIndex_t b);
  void addTriangle(const UnsignedIndex_t a, const UnsignedIndex_t b,
                   const UnsignedIndex_t c);

  PointStorage& getVertexList(void);
  const PointStorage& getVertexList(void) const;

  EdgeStorage& getBoundaryEdgeList(void);
  const EdgeStorage& getBoundaryEdgeList(void) const;

  TriangleStorage& getTriangleList(void);
  const TriangleStorage& getTriangleList(void) const;

  PointStorage::size_type nVertices(void) const;
  EdgeStorage::size_type nBoundaryEdges(void) const;
  TriangleStorage::size_type nTriangles(void) const;

  void clearVertices(void);
  void clearBoundaryEdges(void);
  void clearTriangles(void);
  void clearAll(void);
  void write(const std::string& filename);

  /// Refines triangles to be less than a_max_size. Only two of the dimensions
  /// are used, and the third (a_compute_dim) is comptued according to the
  /// provided function.
  void refineSize(
      const double a_max_size, const UnsignedIndex_t a_compute_dim,
      std::function<double(const double a_x, const double a_y)> a_func);

 private:
  // Includes both edge and hole vertices
  PointStorage vertices_m;
  EdgeStorage bdy_edges_m;
  TriangleStorage triangles_m;
};

inline std::ostream& operator<<(
    std::ostream& out, const TriangulatedSurfaceOutput& a_triangulated_surface);

}  // namespace IRL

#include "irl/surface_mesher/triangulated_surface.tpp"

#endif  // IRL_SURFACE_MESH_TRIANGULATED_SURFACE_H_
