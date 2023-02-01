// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_NORMAL_H_
#define IRL_GEOMETRY_GENERAL_NORMAL_H_

#include <float.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "irl/geometry/general/pt.h"
#include "irl/helpers/byte_buffer.h"
#include "irl/helpers/helper.h"
#include "irl/helpers/mymath.h"
#include "irl/helpers/serializer.h"
#include "irl/parameters/constants.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class ScalarType>
class NormalBase;
/// \brief Overload * operator to be dot product of the two normal vectors
template <class ScalarType>
ScalarType operator*(const NormalBase<ScalarType>& a_normal_0,
                     const NormalBase<ScalarType>& a_normal_1);
/// \brief Overload * operator between PtBase and NormalBase to be a dot
/// product.
template <class ScalarType>
ScalarType operator*(const PtBase<ScalarType>& a_pt,
                     const NormalBase<ScalarType>& a_normal);
/// \brief Overload * operator between PtBase and NormalBase to be a dot
/// product.
template <class ScalarType>
ScalarType operator*(const NormalBase<ScalarType>& a_normal,
                     const PtBase<ScalarType>& a_pt);
/// \brief Overload * operator between double and NormalBase to return a
/// NormalBase.
template <class ScalarType>
NormalBase<ScalarType> operator*(const ScalarType a_double,
                                 const NormalBase<ScalarType>& a_normal);
/// \brief Overload * operator between double and NormalBase to be return a
/// NormalBase.
template <class ScalarType>
NormalBase<ScalarType> operator*(const NormalBase<ScalarType>& a_normal,
                                 const ScalarType a_double);

/// \brief Overload * operator between double and NormalBase to be return a
/// NormalBase.
template <class ScalarType>
NormalBase<ScalarType> operator/(const NormalBase<ScalarType>& a_normal,
                                 const ScalarType a_double);

/// \brief For printing out a normal to a stream.
template <class ScalarType>
inline std::ostream& operator<<(std::ostream& out,
                                const NormalBase<ScalarType>& a_normal);

/// \brief A normal vector in 3D space.
template <class ScalarType>
class NormalBase : public Expr<NormalBase<ScalarType>> {
  using iterator = typename Vec3<ScalarType>::iterator;
  using const_iterator = typename Vec3<ScalarType>::const_iterator;

 public:
  using value_type = ScalarType;
  /// \brief Default constructor with zero normal.
  constexpr NormalBase(void);

  /// \brief Constructor for normal given 3 different values.
  constexpr NormalBase(const UnsignedIndex_t a_normal_x,
                       const UnsignedIndex_t a_normal_y,
                       const UnsignedIndex_t a_normal_z);
  constexpr NormalBase(const int a_normal_x, const int a_normal_y,
                       const int a_normal_z);
  constexpr NormalBase(const double a_normal_x, const double a_normal_y,
                       const double a_normal_z);
  constexpr NormalBase(const Quad_t a_normal_x, const Quad_t a_normal_y,
                       const Quad_t a_normal_z);

  static NormalBase normalized(const ScalarType a_normal_x,
                               const ScalarType a_normal_y,
                               const ScalarType a_normal_z);

  static NormalBase fromRawDoublePointer(const ScalarType* a_normal);

  static NormalBase fromRawDoublePointerNormalized(const ScalarType* a_normal);

  static NormalBase fromPt(const PtBase<ScalarType>& a_pt);

  static NormalBase fromPtNormalized(const PtBase<ScalarType>& a_pt);

  static PtBase<ScalarType> toPt(const NormalBase<ScalarType>& a_normal);

  const NormalBase<double> toDoubleNormal(void) const;
  const NormalBase<Quad_t> toQuadNormal(void) const;

  static NormalBase fromScalarConstant(const ScalarType a_constant);

  // Purposefully not marked explicit to allow elevation from vector_sum to
  // Expr<vector_sum<>>
  template <class E>
  NormalBase(const Expr<E>& a_expr);

  // Purposefully not marked explicit to allow elevation from vector_sum to
  // Expr<vector_sum<>>
  template <class E>
  NormalBase(Expr<E>&& a_expr);

  template <class E>
  NormalBase& operator=(const Expr<E>& a_expr);

  template <class E>
  NormalBase& operator=(Expr<E>&& a_expr);

  /// \brief Return address of `a_d` index of normal.
  ScalarType& operator[](const UnsignedIndex_t a_d);

  /// \brief Const version of `operator[]`.
  const ScalarType& operator[](const UnsignedIndex_t a_d) const;

  /// \brief Length of point vector.
  static constexpr UnsignedIndex_t size(void);

  /// \brief Overload unary - operator.
  NormalBase operator-(void) const;

  /// \brief Increment by another normal.
  NormalBase& operator+=(const NormalBase<ScalarType>& a_rhs);

  /// \brief Overload /= operator to divide all components by scalar.
  NormalBase& operator/=(const double a_value);
  NormalBase& operator/=(const Quad_t a_value);

  /// \brief Overload *= operator to multiply all components by scalar.
  NormalBase& operator*=(const double a_value);
  NormalBase& operator*=(const Quad_t a_value);

  /// \brief Assign all elements in NormalBase the a value.
  NormalBase& operator=(const double a_value);
  NormalBase& operator=(const Quad_t a_value);

  /// \brief Comparison operator to return if two normals are the same
  bool operator==(const NormalBase& a_normal) const;

  /// \brief Comparison operator for const to return if two normals are not the
  /// same.
  bool operator!=(const NormalBase& a_normal) const;

  /// \brief Return the index of the component with maximum magnitude.
  UnsignedIndex_t maxMagnitudeComponent(void) const;

  /// \brief Const version of calculateMagnitude.
  ScalarType calculateMagnitude(void) const;
  ScalarType calculateSquaredMagnitude(void) const;

  /// \brief NormalBaseize the normal
  void normalize(void);

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Return size of the serialized normal class in bytes.
  LargeOffsetIndex_t getSerializedSize(void) const;

  /// \brief Serialize the normal and store in the ByteBuffer
  void serialize(ByteBuffer* a_buffer) const;

  /// \brief Unpack the serialized normal and store.
  void unpackSerialized(ByteBuffer* a_buffer);

  /// \brief Default destructor.
  ~NormalBase(void) = default;

 private:
  /// \brief Constructor for normal given an array containing 3 values.
  explicit constexpr NormalBase(const double* a_normal);
  explicit constexpr NormalBase(const Quad_t* a_normal);

  /// \brief Set a pt as a normal and normalize. Usually occurs
  /// from subtraction of two points.
  explicit NormalBase(const PtBase<ScalarType>& a_pt);

  explicit NormalBase(const double a_constant);
  explicit NormalBase(const Quad_t a_constant);

  /// \brief N_x, N_y, N_z values of normal.
  Vec3<ScalarType> normal_m;
};

using Normal = NormalBase<double>;

}  // namespace IRL

#include "irl/geometry/general/normal.tpp"

#endif  // IRL_GEOMETRY_GENERAL_NORMAL_H_
