// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_PROXY_VERTEX_ACCESS_H_
#define SRC_GEOMETRY_GENERAL_PROXY_VERTEX_ACCESS_H_

#include <algorithm>
#include <initializer_list>

#include "src/geometry/general/pt_list.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class Derived, class GeometryType, UnsignedIndex_t kNumberOfVertices>
class ProxyVertexAccess {
 public:
  static constexpr UnsignedIndex_t number_of_vertices = kNumberOfVertices;
  using VertexType = typename GeometryType::pt_type;

  ProxyVertexAccess(const GeometryType& a_base_geometry,
                    std::initializer_list<UnsignedIndex_t> a_list);

  ProxyVertexAccess(const GeometryType& a_base_geometry,
                    std::array<UnsignedIndex_t, kNumberOfVertices> a_list);

  const VertexType& access(const UnsignedIndex_t a_index) const;

  const std::array<UnsignedIndex_t, kNumberOfVertices>& getIndexMapping(
      void) const;

  /// \brief Return the number of vertices in this polygon.
  static constexpr UnsignedIndex_t getNumberOfVerticesInObject(void);

  /// \brief Return size of the serialized Polyhedron.
  LargeOffsetIndex_t getSerializedSize(void) const;

  /// \brief Serialize and pack the Polyhedron, copying the points.
  /// This will need to be unpacked into an object of this
  /// type that has its own storage.
  void serialize(ByteBuffer* a_buffer) const;

  ProxyVertexAccess(const ProxyVertexAccess& a_other);
  ProxyVertexAccess& operator=(const ProxyVertexAccess& a_other);
  ProxyVertexAccess(ProxyVertexAccess&& a_other);
  ProxyVertexAccess& operator=(ProxyVertexAccess&& a_other);

 private:
  const GeometryType& base_geometry_m;
  std::array<UnsignedIndex_t, kNumberOfVertices> index_mapping_m;
};

}  // namespace IRL

#include "src/geometry/general/proxy_vertex_access.tpp"

#endif  // SRC_GEOMETRY_GENERAL_PROXY_VERTEX_ACCESS_H_
