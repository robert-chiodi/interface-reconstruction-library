// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2023 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_PT_H_
#define IRL_GEOMETRY_GENERAL_PT_H_

#include <quadmath.h>
#include <array>
#include <cassert>
#include <iomanip>
#include <ostream>
#include <type_traits>

#include "irl/geometry/general/math_vector.h"
#include "irl/helpers/byte_buffer.h"
#include "irl/helpers/expression_templates.h"
#include "irl/helpers/mymath.h"
#include "irl/helpers/serializer.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class ScalarType>
class PtBase;

template <class ScalarType>
inline std::ostream& operator<<(std::ostream& out,
                                const PtBase<ScalarType>& a_pt);

/// \brief A point in 3D space.
template <class ScalarType>
class PtBase : public Expr<PtBase<ScalarType>> {
 public:
  using value_type = ScalarType;
  /// \brief The default constructor, performs NO INITIALIZATION.
  constexpr PtBase(void);

  /// \brief A constructor that initializes the point given the x,y, and z
  /// locations.
  constexpr PtBase(const UnsignedIndex_t a_x, const UnsignedIndex_t a_y,
                   const UnsignedIndex_t a_z);
  constexpr PtBase(const int a_x, const int a_y, const int a_z);
  constexpr PtBase(const double a_x, const double a_y, const double a_z);
  constexpr PtBase(const Quad_t a_x, const Quad_t a_y, const Quad_t a_z);

  /// \brief Performs the interpolation on an edge between points with known
  /// distances to the plane to calculate the intersection point.
  ///
  /// This function takes two vertices of an edge and the corresponding signed
  /// distances to an intersecting plane. The point of intersection (where
  /// dist=0) is then calculated and returned.
  ///
  /// \param[in] a_point_1 First point on the edge.
  /// \param[in] a_dist_1 Corresponding signed distance from intersecting
  ///  plane for a_point_1.
  /// \param[in] a_point_2 Second point on the edge.
  /// \param[in] a_dist_2 Corresponding signed distance from intersecting
  ///  plane for a_point_2.
  static PtBase fromEdgeIntersection(const PtBase& a_pt_0,
                                     const ScalarType a_dist_0,
                                     const PtBase& a_pt_1,
                                     const ScalarType a_dist_1);

  static PtBase fromRawDoublePointer(const ScalarType* a_loc);

  static PtBase fromArray(const std::array<ScalarType, 3>& a_loc);

  static constexpr PtBase fromScalarConstant(const ScalarType a_value);

  static PtBase fromVec3(const Vec3<ScalarType>& a_vec);

  operator Vec3<ScalarType>(void);

  /// \brief Assignment of point to a constant.
  PtBase& operator=(const ScalarType a_value);

  /// \brief Multiply point by constant value.
  PtBase& operator*=(const ScalarType a_value);

  /// \brief Functions necessary to have consistent interface to PtWithData
  PtBase& getPt(void);
  const PtBase& getPt(void) const;

  /// \brief Functions necessary to have convert to other precisions
  const PtBase<double> toDoublePt(void) const;
  const PtBase<Quad_t> toQuadPt(void) const;

  /// \brief Length of point vector.
  template <class C>
  friend constexpr UnsignedIndex_t size(const PtBase<C>& x);

  /// \brief Return reference to `a_d` index of the point.
  ScalarType& operator[](const UnsignedIndex_t a_d);

  /// \brief Const version of `operator[]`.
  const ScalarType& operator[](const UnsignedIndex_t a_d) const;

  /// \brief Return value of `loc_m[0]` (x location).
  ScalarType& x(void);
  /// \brief Return value of `loc_m[1]` (y location).
  ScalarType& y(void);
  /// \brief Return value of `loc_m[2]` (z location).
  ScalarType& z(void);
  /// \brief Const version of `x()`.
  const ScalarType& x(void) const;
  /// \brief Const version of `y()`.
  const ScalarType& y(void) const;
  /// \brief Const version of `z()`.
  const ScalarType& z(void) const;

  /// \brief Return the index of the dimension with maximum magnitude.
  UnsignedIndex_t maxComponent(void) const;

  // \brief Unary minus operator
  PtBase operator-(void) const;

  /// \brief Increment this point by another point.
  PtBase& operator+=(const PtBase<ScalarType>& a_rhs);

  /// \brief Decrement this point by another point.
  PtBase& operator-=(const PtBase<ScalarType>& a_rhs);

  /// \brief Divide each point location by the constant `a_rhs`.
  PtBase operator/(const ScalarType a_rhs);

  /// \brief Divide each point location by the constant `a_rhs`.
  PtBase& operator/=(const ScalarType a_rhs);

  /// \brief Return size of the serialized point class in bytes.
  LargeOffsetIndex_t getSerializedSize(void) const;

  /// \brief Serialize the point and store in the ByteBuffer
  void serialize(ByteBuffer* a_buffer) const;

  /// \brief Unpack the serialized point and store.
  void unpackSerialized(ByteBuffer* a_buffer);

  // Purposefully not marked explicit to allow elevation from vector_sum to
  // Expr<vector_sum<>>
  template <class E>
  PtBase(const Expr<E>& a_expr);

  // Purposefully not marked explicit to allow elevation from vector_sum to
  // Expr<vector_sum<>>
  template <class E>
  PtBase(Expr<E>&& a_expr);

  template <class E>
  PtBase& operator=(const Expr<E>& a_expr);

  template <class E>
  PtBase& operator=(Expr<E>&& a_expr);

  /// \brief Default destructor.
  ~PtBase(void) = default;

 private:
  /// \brief A constructor that initializes the
  /// point given 3 elements in the array `a_loc`.
  explicit constexpr PtBase(const ScalarType* a_loc);

  /// \brief A constructor that initializes the
  /// point given 3 elements in the array `a_loc`.
  explicit constexpr PtBase(const std::array<ScalarType, 3>& a_loc);

  /// \brief Construct from a scalar constant value.
  explicit constexpr PtBase(const ScalarType a_value);

  /// \brief Construct from a Vec3 MathVector.
  explicit PtBase(const Vec3<ScalarType>& a_vec);

  std::array<ScalarType, 3> loc_m;  ///< x,y,z (0,1,2) location of the point.
};

using Pt = PtBase<double>;

}  // namespace IRL

#include "irl/geometry/general/pt.tpp"

#endif  // IRL_GEOMETRY_GENERAL_PT_H_
