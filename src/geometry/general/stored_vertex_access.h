// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_STORED_VERTEX_ACCESS_H_
#define SRC_GEOMETRY_GENERAL_STORED_VERTEX_ACCESS_H_

#include <initializer_list>

#include "src/geometry/general/pt_list.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class Derived, class VertexType, UnsignedIndex_t kNumberOfVertices>
class StoredVertexAccess {
 public:
  static constexpr UnsignedIndex_t number_of_vertices = kNumberOfVertices;

  StoredVertexAccess(void) = default;

  StoredVertexAccess(std::initializer_list<VertexType> a_list);

  /// \brief Construct n-pts form array of points.
  static Derived fromRawPtPointer(const UnsignedIndex_t a_number_of_pts,
                                  const VertexType* a_array_of_pts);

  static Derived fromRawDoublePointer(const UnsignedIndex_t a_number_of_pts,
                                      const double* a_array_of_locs);

  VertexType& access(const UnsignedIndex_t a_index);

  const VertexType& access(const UnsignedIndex_t a_index) const;

  /// \brief Return the number of vertices in this polygon.
  static constexpr UnsignedIndex_t getNumberOfVerticesInObject(void);

  /// \brief Return size of the serialized Polyhedron.
  LargeOffsetIndex_t getSerializedSize(void) const;

  /// \brief Serialize and pack the Polyhedron.
  void serialize(ByteBuffer* a_buffer) const;

  /// \brief Unpack the polyhedron.
  void unpackSerialized(ByteBuffer* a_buffer);

 protected:
  /// \brief Construct n-pts form array of pts.
  StoredVertexAccess(const UnsignedIndex_t a_number_of_pts,
                     const VertexType* a_array_of_pts);

  /// \brief Construct n-pts form array of doubles.
  StoredVertexAccess(const UnsignedIndex_t a_number_of_pts,
                     const double* a_array_of_locs);

 private:
  PtList<VertexType, kNumberOfVertices> pt_list_m;
};

}  // namespace IRL

#include "src/geometry/general/stored_vertex_access.tpp"

#endif  // SRC_GEOMETRY_GENERAL_STORED_VERTEX_ACCESS_H_
