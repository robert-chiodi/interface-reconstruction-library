// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_MATH_VECTOR_H_
#define SRC_GEOMETRY_GENERAL_MATH_VECTOR_H_

#include <algorithm>
#include <cassert>

#include "src/helpers/byte_buffer.h"
#include "src/helpers/mymath.h"
#include "src/helpers/serializer.h"
#include "src/parameters/defined_types.h"
#include "src/helpers/expression_templates.h"

namespace IRL {

template <UnsignedIndex_t kNumberOfElements>
class MathVector : public Expr<MathVector<kNumberOfElements>>{
 public:
  using value_type = double;
  using iterator = typename std::array<double, kNumberOfElements>::iterator;
  using const_iterator =
      typename std::array<double, kNumberOfElements>::const_iterator;

  /// \brief Default constructor.
  constexpr MathVector(void) = default;

  /// \brief Constructor for vector given 3 different values.
  constexpr MathVector(const double a_element_0, const double a_element_1,
                       const double a_element_2);
  static MathVector fromRawDoublePointer(const double* a_vec);

  static MathVector fromScalarConstant(const double a_constant);

  // Purposefully not marked explicit to allow elevation from vector_sum to
  // Expr<vector_sum<>>
  template <class E>
  MathVector(const Expr<E>& a_expr);

  // Purposefully not marked explicit to allow elevation from vector_sum to
  // Expr<vector_sum<>>
  template <class E>
  MathVector(Expr<E>&& a_expr);

  template <class E>
  MathVector& operator=(const Expr<E>& a_expr);

  template <class E>
  MathVector& operator=(Expr<E>&& a_expr);

  /// \brief Length of point vector.
  static constexpr UnsignedIndex_t size(void);

  /// \brief Return address of `a_d` index of the MathVector.
  double& operator[](const UnsignedIndex_t a_d);

  /// \brief Const version of `operator[]`.
  const double& operator[](const UnsignedIndex_t a_d) const;

  MathVector operator-(void) const;

  MathVector& operator+=(const MathVector& a_other_vector);

  MathVector& operator-=(const MathVector& a_other_vector);

  MathVector& operator+=(const double a_scalar);

  MathVector& operator-=(const double a_scalar);

  MathVector& operator*=(const double a_scalar);

  MathVector& operator/=(const double a_scalar);

  /// \brief Return the index of the component with maximum magnitude.
  UnsignedIndex_t calculateIndexOfLargestMagnitude(void) const;

  double calculateMagnitude(void) const;

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

 private:
  /// \brief Constructor for normal given an array containing 3 values.
  explicit MathVector(const double* a_vec);

  explicit MathVector(const double a_constant);

  std::array<double, kNumberOfElements> elements_m;
};

using Vec3 = MathVector<3>;

/// \brief Overload * operator to be dot product of the two vectors
template <UnsignedIndex_t kNumberOfElements>
__attribute__((const)) inline double operator*(
    const MathVector<kNumberOfElements>& a_vec_0,
    const MathVector<kNumberOfElements>& a_vec_1);

/// \brief Overload + operator between MathVectors of same length to add
/// elements.
template <UnsignedIndex_t kNumberOfElements>
__attribute__((const)) inline MathVector<kNumberOfElements> operator+(
    const MathVector<kNumberOfElements>& a_vec_0,
    const MathVector<kNumberOfElements>& a_vec_1);

/// \brief Overload - operator between MathVectors of same length to add
/// elements.
template <UnsignedIndex_t kNumberOfElements>
__attribute__((const)) inline MathVector<kNumberOfElements> operator-(
    const MathVector<kNumberOfElements>& a_vec_0,
    const MathVector<kNumberOfElements>& a_vec_1);

/// \brief Overload * operator between double and MathVector to return a
/// MathVector.
template <UnsignedIndex_t kNumberOfElements>
__attribute__((const)) inline MathVector<kNumberOfElements> operator*(
    const double a_double, const MathVector<kNumberOfElements>& a_vec);
/// \brief Overload * operator between double and MathVector to be return a
/// MathVector.
template <UnsignedIndex_t kNumberOfElements>
__attribute__((const)) inline MathVector<kNumberOfElements> operator*(
    const MathVector<kNumberOfElements>& a_vec, const double a_double);

/// \brief Overload / operator between double and MathVector to be return a
/// MathVector.
template <UnsignedIndex_t kNumberOfElements>
__attribute__((const)) inline MathVector<kNumberOfElements> operator/(
    const MathVector<kNumberOfElements>& a_vec, const double a_double);

}  // namespace IRL

#include "src/geometry/general/math_vector.tpp"

#endif  // SRC_GEOMETRY_GENERAL_MATH_VECTOR_H_
