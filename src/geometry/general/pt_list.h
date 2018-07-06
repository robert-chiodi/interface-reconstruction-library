// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_PT_LIST_H_
#define SRC_GEOMETRY_GENERAL_PT_LIST_H_

#include <algorithm>
#include <array>
#include <cstring>
#include <initializer_list>

#include "src/geometry/general/pt.h"
#include "src/helpers/byte_buffer.h"
#include "src/helpers/serializer.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
class Polyhedron;

template <class VertexType, UnsignedIndex_t kMaxNumberOfPts>
class PtList {
  friend class Polyhedron<VertexType, kMaxNumberOfPts>;

 public:
  using iterator = typename std::array<VertexType, kMaxNumberOfPts>::iterator;
  using const_iterator =
      typename std::array<VertexType, kMaxNumberOfPts>::const_iterator;

  /// \brief Default constructor
  PtList(void) = default;

  PtList(std::initializer_list<VertexType> a_list);

  static PtList fromRawPtPointer(const UnsignedIndex_t a_number_of_pts,
                                 const VertexType* a_array_of_pts);

  static PtList fromRawDoublePointer(const UnsignedIndex_t a_number_of_pts,
                                     const double* a_array_of_locs);

  /// \brief Return the number of vertices in this polygon.
  static constexpr UnsignedIndex_t getNumberOfPts(void);

  /// \brief Return const pointer to the vertices
  const VertexType* getPtList(void) const;

  /// \brief Access through overloaded operator[]
  VertexType& operator[](const UnsignedIndex_t a_index);

  /// \brief Const access through overloaded operator[]
  const VertexType& operator[](const UnsignedIndex_t a_index) const;

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

  /// \brief Return size of the serialized Ptist.
  LargeOffsetIndex_t getSerializedSize(void) const;

  /// \brief Serialize and pack the points.
  void serialize(ByteBuffer* a_buffer) const;

  /// \brief Unpack the points and store.
  void unpackSerialized(ByteBuffer* a_buffer);
  /// \brief Default destructor
  ~PtList(void) = default;

 protected:
  /// \brief Construct a tetrahedron
  constexpr PtList(const VertexType& a_pt0, const VertexType& a_pt1,
                   const VertexType& a_pt2, const VertexType& a_pt3);
  /// \brief Construct n-pts form array of points.
  PtList(const UnsignedIndex_t a_number_of_pts,
         const VertexType* a_array_of_pts);

  /// \brief Construct n-pts form array of doubles.
  PtList(const UnsignedIndex_t a_number_of_pts, const double* a_array_of_locs);

 private:
  std::array<VertexType, kMaxNumberOfPts> pt_list_m;
};

}  // namespace IRL

#include "src/geometry/general/pt_list.tpp"

#endif  // SRC_GEOMETRY_GENERAL_PT_LIST_H_
