// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_EXPANDABLE_PT_LIST_H_
#define IRL_GEOMETRY_GENERAL_EXPANDABLE_PT_LIST_H_

#include <algorithm>
#include <cstring>
#include <iostream>

#include "irl/data_structures/small_vector.h"
#include "irl/geometry/general/pt.h"
#include "irl/helpers/byte_buffer.h"
#include "irl/helpers/serializer.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class VertexType>
class PolygonBase;

template <class VertexType, UnsignedIndex_t kStaticAllocSize = 10>
class ExpandablePtList {
  friend class PolygonBase<VertexType>;

 public:
  static constexpr UnsignedIndex_t on_stack_size = kStaticAllocSize;
  using StorageType = SmallVector<VertexType, on_stack_size>;
  using iterator = typename StorageType::iterator;
  using const_iterator = typename StorageType::const_iterator;
  /// \brief Default constructor
  ExpandablePtList(void) = default;

  static ExpandablePtList fromRawPtPointer(
      const UnsignedIndex_t a_number_of_pts, const VertexType* a_array_of_pts);

  static ExpandablePtList fromRawDoublePointer(
      const UnsignedIndex_t a_number_of_pts, const double* a_array_of_locs);
  /// \brief Const return the number of vertices in this polygon.
  UnsignedIndex_t getNumberOfPts(void) const;

  /// \brief Set the number of vertices directly.
  void setNumberOfPts(const UnsignedIndex_t a_number);

  /// \brief Reserve space in the Pt vector.
  void reserve(const UnsignedIndex_t a_number);

  /// \brief Add point to the point list.
  void addPt(const VertexType& a_pt);

  /// \brief Add point to the point list.
  void addPtAtIndex(const UnsignedIndex_t a_index, const VertexType& a_pt);

  /// \brief Access through overloaded operator[]
  VertexType& operator[](const UnsignedIndex_t a_index);

  /// \brief Const access through overloaded operator[]
  const VertexType& operator[](const UnsignedIndex_t a_index) const;

  /// \brief Remove point from the point list.
  void removePt(const UnsignedIndex_t a_index);

  /// \brief Remove last point from the point list.
  void removeLastPt(void);

  /// \brief Return a point for the lower limits of the polygon in 3D space.
  IRL::Pt getLowerLimits(void) const;

  /// \brief Return a point for the upper limits of the polygon in 3D space.
  IRL::Pt getUpperLimits(void) const;

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Return size of the serialized ExpandablePtList.
  LargeOffsetIndex_t getSerializedSize(void) const;

  /// \brief Serialize and pack the points.
  void serialize(ByteBuffer* a_buffer) const;

  /// \brief Unpack the points and store
  void unpackSerialized(ByteBuffer* a_buffer);
  /// \brief Default destructor
  ~ExpandablePtList(void) = default;

 protected:
  /// \brief Construct n-pts form array of points.
  ExpandablePtList(const UnsignedIndex_t a_number_of_pts,
                   const VertexType* a_array_of_pts);

  /// \brief Construct n-pts form array of doubles.
  ExpandablePtList(const UnsignedIndex_t a_number_of_pts,
                   const double* a_array_of_locs);

 private:
  void checkIfStaticAllocationExceeded(void) const;
  StorageType pt_list_m;
};
}  // namespace IRL

#include "irl/geometry/general/expandable_pt_list.tpp"

#endif // IRL_GEOMETRY_GENERAL_EXPANDABLE_PT_LIST_H_
