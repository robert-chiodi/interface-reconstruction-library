// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_DELAYED_EXPANDABLE_VERTEX_ACCESS_H_
#define SRC_GEOMETRY_GENERAL_DELAYED_EXPANDABLE_VERTEX_ACCESS_H_

#include <initializer_list>

#include "src/geometry/general/expandable_pt_list.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class Derived, class VertexType, UnsignedIndex_t kMaskedVertices>
class DelayedExpandableVertexAccess {
 public:
  DelayedExpandableVertexAccess(void);

  DelayedExpandableVertexAccess(std::initializer_list<VertexType> a_list);

  /// \brief Construct n-pts form array of points.
  static Derived fromRawPtPointer(const UnsignedIndex_t a_number_of_pts,
                                  const VertexType* a_array_of_pts);

  static Derived fromRawDoublePointer(const UnsignedIndex_t a_number_of_pts,
                                      const double* a_array_of_locs);

  VertexType& access(const UnsignedIndex_t a_index);

  const VertexType& access(const UnsignedIndex_t a_index) const;

  VertexType& unmaskedAccess(const UnsignedIndex_t a_index);

  const VertexType& unmaskedAccess(const UnsignedIndex_t a_index) const;

  /// \brief Return the number of vertices in this polygon.
  UnsignedIndex_t getNumberOfVerticesInObject(void) const;

  void setNumberOfVerticesInObject(const UnsignedIndex_t a_size);

  void addVertex(const VertexType& a_vertex);
  void addVertexAtIndex(const UnsignedIndex_t a_index, const VertexType& a_pt);

  void reserve(const UnsignedIndex_t a_size);

  /// \brief Add point to the point list.
  void removeLastVertex(void);

  void removePt(const UnsignedIndex_t a_index);

  /// \brief Return size of the serialized Polyhedron.
  LargeOffsetIndex_t getSerializedSize(void) const;

  /// \brief Serialize and pack the Polyhedron.
  void serialize(ByteBuffer* a_buffer) const;

  /// \brief Unpack the polyhedron.
  void unpackSerialized(ByteBuffer* a_buffer);

 protected:
  /// \brief Construct n-pts form array of pts.
  DelayedExpandableVertexAccess(const UnsignedIndex_t a_number_of_pts,
                                const VertexType* a_array_of_pts);

  /// \brief Construct n-pts form array of doubles.
  DelayedExpandableVertexAccess(const UnsignedIndex_t a_number_of_pts,
                                const double* a_array_of_locs);
  void insertMaskedVertices(void);

 private:
  ExpandablePtList<VertexType> pt_list_m;
};

// Used to remove masking by shadowing [] operator to give maskedAccess instead
// of access.
template <class MaskedType>
class MaskStripper : public MaskedType {
 public:
  using pt_type = typename MaskedType::pt_type;

  pt_type& operator[](const UnsignedIndex_t a_index);
  const pt_type& operator[](const UnsignedIndex_t a_index) const;
};

}  // namespace IRL

#include "src/geometry/general/delayed_expandable_vertex_access.tpp"

#endif  // SRC_GEOMETRY_GENERAL_DELAYED_EXPANDABLE_VERTEX_ACCESS_H_
