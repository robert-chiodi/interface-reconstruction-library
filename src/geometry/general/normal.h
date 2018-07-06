// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_NORMAL_H_
#define SRC_GEOMETRY_GENERAL_NORMAL_H_

#include <float.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "src/geometry/general/pt.h"
#include "src/helpers/byte_buffer.h"
#include "src/helpers/helper.h"
#include "src/helpers/mymath.h"
#include "src/helpers/serializer.h"
#include "src/parameters/constants.h"
#include "src/parameters/defined_types.h"

namespace IRL {

class Normal;
/// \brief Overload * operator to be dot product of the two normal vectors
double operator*(const Normal& a_normal_0, const Normal& a_normal_1);
/// \brief Overload * operator between Pt and Normal to be a dot product.
double operator*(const Pt& a_pt, const Normal& a_normal);
/// \brief Overload * operator between Pt and Normal to be a dot product.
double operator*(const Normal& a_normal, const Pt& a_pt);
/// \brief Overload * operator between double and Normal to return a Normal.
Normal operator*(const double a_double, const Normal& a_normal);
/// \brief Overload * operator between double and Normal to be return a Normal.
Normal operator*(const Normal& a_normal, const double a_double);

/// \brief Overload * operator between double and Normal to be return a Normal.
Normal operator/(const Normal& a_normal, const double a_double);

/// \brief For printing out a normal to a stream.
inline std::ostream& operator<<(std::ostream& out, const Normal& a_normal);

/// \brief A normal vector in 3D space.
class Normal : public Expr<Normal>{
  using iterator = Vec3::iterator;
  using const_iterator = Vec3::const_iterator;

 public:
  using value_type = double;
  /// \brief Default constructor with zero normal.
  constexpr Normal(void);

  /// \brief Constructor for normal given 3 different values.
  constexpr Normal(const double a_normal_x, const double a_normal_y,
                   const double a_normal_z);

  static Normal normalized(const double a_normal_x, const double a_normal_y,
                           const double a_normal_z);

  static Normal fromRawDoublePointer(const double* a_normal);

  static Normal fromRawDoublePointerNormalized(const double* a_normal);

  static Normal fromPt(const Pt& a_pt);

  static Normal fromPtNormalized(const Pt& a_pt);

  static Pt toPt(const Normal& a_normal);

  static Normal fromScalarConstant(const double a_constant);

  // Purposefully not marked explicit to allow elevation from vector_sum to
  // Expr<vector_sum<>>
  template <class E>
  Normal(const Expr<E>& a_expr);

  // Purposefully not marked explicit to allow elevation from vector_sum to
  // Expr<vector_sum<>>
  template <class E>
  Normal(Expr<E>&& a_expr);

  template <class E>
  Normal& operator=(const Expr<E>& a_expr);

  template <class E>
  Normal& operator=(Expr<E>&& a_expr);

  /// \brief Return address of `a_d` index of normal.
  double& operator[](const UnsignedIndex_t a_d);

  /// \brief Const version of `operator[]`.
  const double& operator[](const UnsignedIndex_t a_d) const;

  /// \brief Length of point vector.
  static constexpr UnsignedIndex_t size(void);

  /// \brief Overload unary - operator.
  Normal operator-(void) const;

  /// \brief Increment by another normal.
  Normal& operator+=(const Normal& a_rhs);

  /// \brief Overload /= operator to divide all components by scalar.
  Normal& operator/=(const double a_value);

  /// \brief Overload *= operator to multiply all components by scalar.
  Normal& operator*=(const double a_value);

  /// \brief Assign all elements in Normal the a value.
  Normal& operator=(const double a_value);

  /// \brief Comparison operator to return if two normals are the same
  bool operator==(const Normal& a_normal) const;

  /// \brief Comparison operator for const to return if two normals are not the
  /// same.
  bool operator!=(const Normal& a_normal) const;

  /// \brief Return the index of the component with maximum magnitude.
  UnsignedIndex_t maxMagnitudeComponent(void) const;

  /// \brief Const version of calculateMagnitude.
  double calculateMagnitude(void) const;

  /// \brief Normalize the normal
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
  ~Normal(void) = default;

 private:
  /// \brief Constructor for normal given an array containing 3 values.
  explicit constexpr Normal(const double* a_normal);

  /// \brief Set a pt as a normal and normalize. Usually occurs
  /// from subtraction of two points.
  explicit Normal(const Pt& a_pt);

  explicit Normal(const double a_constant);

  /// \brief N_x, N_y, N_z values of normal.
  Vec3 normal_m;
};

}  // namespace IRL

#include "src/geometry/general/normal.tpp"

#endif  // SRC_GEOMETRY_GENERAL_NORMAL_H_
