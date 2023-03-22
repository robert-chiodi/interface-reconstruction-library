// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_MATH_VECTOR_H_
#define IRL_GEOMETRY_GENERAL_MATH_VECTOR_H_

#include <algorithm>
#include <cassert>

#include "irl/helpers/byte_buffer.h"
#include "irl/helpers/mymath.h"
#include "irl/helpers/serializer.h"
#include "irl/parameters/defined_types.h"

#include "irl/helpers/expression_templates.h"

namespace IRL {

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
class MathVector : public Expr<MathVector<kNumberOfElements, ScalarType>> {
 public:
  using value_type = ScalarType;
  using iterator = typename std::array<ScalarType, kNumberOfElements>::iterator;
  using const_iterator =
      typename std::array<ScalarType, kNumberOfElements>::const_iterator;

  /// \brief Default constructor.
  constexpr MathVector(void) = default;

  /// \brief Constructor for vector given 3 different values.
  constexpr MathVector(const ScalarType a_element_0,
                       const ScalarType a_element_1,
                       const ScalarType a_element_2);
  static MathVector fromRawDoublePointer(const ScalarType* a_vec);

  static MathVector fromScalarConstant(const ScalarType a_constant);

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
  ScalarType& operator[](const UnsignedIndex_t a_d);

  /// \brief Const version of `operator[]`.
  const ScalarType& operator[](const UnsignedIndex_t a_d) const;

  MathVector operator-(void) const;

  MathVector& operator+=(const MathVector& a_other_vector);

  MathVector& operator-=(const MathVector& a_other_vector);

  MathVector& operator+=(const ScalarType a_scalar);

  MathVector& operator-=(const ScalarType a_scalar);

  MathVector& operator*=(const ScalarType a_scalar);

  MathVector& operator/=(const ScalarType a_scalar);

  /// \brief Return the index of the component with maximum magnitude.
  UnsignedIndex_t calculateIndexOfLargestMagnitude(void) const;

  ScalarType calculateMagnitude(void) const;

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
  explicit MathVector(const ScalarType* a_vec);

  explicit MathVector(const ScalarType a_constant);

  std::array<ScalarType, kNumberOfElements> elements_m;
};

template <class ScalarType>
using Vec3 = MathVector<3, ScalarType>;

/// \brief Overload * operator to be dot product of the two vectors
template <UnsignedIndex_t kNumberOfElements, class ScalarType>
__attribute__((const)) inline ScalarType operator*(
    const MathVector<kNumberOfElements, ScalarType>& a_vec_0,
    const MathVector<kNumberOfElements, ScalarType>& a_vec_1);

/// \brief Overload + operator between MathVectors of same length to add
/// elements.
template <UnsignedIndex_t kNumberOfElements, class ScalarType>
__attribute__((const)) inline MathVector<kNumberOfElements, ScalarType>
operator+(const MathVector<kNumberOfElements, ScalarType>& a_vec_0,
          const MathVector<kNumberOfElements, ScalarType>& a_vec_1);

/// \brief Overload - operator between MathVectors of same length to add
/// elements.
template <UnsignedIndex_t kNumberOfElements, class ScalarType>
__attribute__((const)) inline MathVector<kNumberOfElements, ScalarType>
operator-(const MathVector<kNumberOfElements, ScalarType>& a_vec_0,
          const MathVector<kNumberOfElements, ScalarType>& a_vec_1);

/// \brief Overload * operator between scalar and MathVector to return a
/// MathVector.
template <UnsignedIndex_t kNumberOfElements, class ScalarType>
__attribute__((const)) inline MathVector<kNumberOfElements, ScalarType>
operator*(const ScalarType a_scalar,
          const MathVector<kNumberOfElements, ScalarType>& a_vec);
/// \brief Overload * operator between scalar and MathVector to be return a
/// MathVector.
template <UnsignedIndex_t kNumberOfElements, class ScalarType>
__attribute__((const)) inline MathVector<kNumberOfElements, ScalarType>
operator*(const MathVector<kNumberOfElements, ScalarType>& a_vec,
          const ScalarType a_scalar);

/// \brief Overload / operator between scalar and MathVector to be return a
/// MathVector.
template <UnsignedIndex_t kNumberOfElements, class ScalarType>
__attribute__((const)) inline MathVector<kNumberOfElements, ScalarType>
operator/(const MathVector<kNumberOfElements, ScalarType>& a_vec,
          const ScalarType a_scalar);

}  // namespace IRL

#include "irl/geometry/general/math_vector.tpp"

#endif  // IRL_GEOMETRY_GENERAL_MATH_VECTOR_H_
