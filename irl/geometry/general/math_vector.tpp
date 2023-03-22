// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_MATH_VECTOR_TPP_
#define IRL_GEOMETRY_GENERAL_MATH_VECTOR_TPP_

#include <algorithm>

namespace IRL {

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
constexpr MathVector<kNumberOfElements, ScalarType>::MathVector(
    const ScalarType a_element_0, const ScalarType a_element_1,
    const ScalarType a_element_2)
    : elements_m{a_element_0, a_element_1, a_element_2} {}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType>
MathVector<kNumberOfElements, ScalarType>::fromRawDoublePointer(
    const ScalarType* a_vec) {
  return MathVector(a_vec);
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType>
MathVector<kNumberOfElements, ScalarType>::fromScalarConstant(
    const ScalarType a_constant) {
  return MathVector(a_constant);
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
template <class E>
inline MathVector<kNumberOfElements, ScalarType>::MathVector(
    const Expr<E>& a_expr) {
  const E& expr(a_expr);
  for (UnsignedIndex_t n = 0; n < kNumberOfElements; ++n) {
    elements_m[n] = expr[n];
  }
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
template <class E>
inline MathVector<kNumberOfElements, ScalarType>::MathVector(Expr<E>&& a_expr) {
  const E& expr(a_expr);
  for (UnsignedIndex_t n = 0; n < kNumberOfElements; ++n) {
    elements_m[n] = expr[n];
  }
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
template <class E>
inline MathVector<kNumberOfElements, ScalarType>&
MathVector<kNumberOfElements, ScalarType>::operator=(const Expr<E>& a_expr) {
  const E& expr(a_expr);
  for (UnsignedIndex_t n = 0; n < kNumberOfElements; ++n) {
    elements_m[n] = expr[n];
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
template <class E>
inline MathVector<kNumberOfElements, ScalarType>&
MathVector<kNumberOfElements, ScalarType>::operator=(Expr<E>&& a_expr) {
  const E& expr(a_expr);
  for (UnsignedIndex_t n = 0; n < kNumberOfElements; ++n) {
    elements_m[n] = expr[n];
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
constexpr UnsignedIndex_t MathVector<kNumberOfElements, ScalarType>::size(
    void) {
  return kNumberOfElements;
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
ScalarType& MathVector<kNumberOfElements, ScalarType>::operator[](
    const UnsignedIndex_t a_d) {
  assert(a_d < this->size());
  return elements_m[a_d];
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
const ScalarType& MathVector<kNumberOfElements, ScalarType>::operator[](
    const UnsignedIndex_t a_d) const {
  assert(a_d < this->size());
  return elements_m[a_d];
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType>
MathVector<kNumberOfElements, ScalarType>::operator-(void) const {
  MathVector vector_to_return = (*this);
  for (auto& element : (*this)) {
    element = -element;
  }
  return vector_to_return;
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType>&
MathVector<kNumberOfElements, ScalarType>::operator+=(
    const MathVector& a_other_vector) {
  assert(a_other_vector.size() == this->size());
  for (UnsignedIndex_t ind = 0; ind < a_other_vector.size(); ++ind) {
    (*this)[ind] += a_other_vector[ind];
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType>&
MathVector<kNumberOfElements, ScalarType>::operator-=(
    const MathVector& a_other_vector) {
  assert(a_other_vector.size() == this->size());
  for (UnsignedIndex_t ind = 0; ind < a_other_vector.size(); ++ind) {
    (*this)[ind] -= a_other_vector[ind];
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType>&
MathVector<kNumberOfElements, ScalarType>::operator+=(
    const ScalarType a_scalar) {
  for (auto& element : (*this)) {
    element += a_scalar;
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType>&
MathVector<kNumberOfElements, ScalarType>::operator-=(
    const ScalarType a_scalar) {
  for (auto& element : (*this)) {
    element -= a_scalar;
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType>&
MathVector<kNumberOfElements, ScalarType>::operator*=(
    const ScalarType a_scalar) {
  for (auto& element : (*this)) {
    element *= a_scalar;
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType>&
MathVector<kNumberOfElements, ScalarType>::operator/=(
    const ScalarType a_scalar) {
  for (auto& element : (*this)) {
    element /= a_scalar;
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
UnsignedIndex_t
MathVector<kNumberOfElements, ScalarType>::calculateIndexOfLargestMagnitude(
    void) const {
  UnsignedIndex_t index_of_max_magnitude;
  ScalarType max_value = -DBL_MAX;
  ScalarType abs_value_of_element;
  for (UnsignedIndex_t ind = 0; ind < this->size(); ++ind) {
    abs_value_of_element = fabs(elements_m[ind]);
    if (abs_value_of_element > max_value) {
      max_value = abs_value_of_element;
      index_of_max_magnitude = ind;
    }
  }
  return index_of_max_magnitude;
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
ScalarType MathVector<kNumberOfElements, ScalarType>::calculateMagnitude(
    void) const {
  return magnitude((*this));
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
typename MathVector<kNumberOfElements, ScalarType>::iterator
MathVector<kNumberOfElements, ScalarType>::begin(void) noexcept {
  return elements_m.begin();
}
template <UnsignedIndex_t kNumberOfElements, class ScalarType>
typename MathVector<kNumberOfElements, ScalarType>::const_iterator
MathVector<kNumberOfElements, ScalarType>::begin(void) const noexcept {
  return this->cbegin();
}
template <UnsignedIndex_t kNumberOfElements, class ScalarType>
typename MathVector<kNumberOfElements, ScalarType>::const_iterator
MathVector<kNumberOfElements, ScalarType>::cbegin(void) const noexcept {
  return elements_m.cbegin();
}
template <UnsignedIndex_t kNumberOfElements, class ScalarType>
typename MathVector<kNumberOfElements, ScalarType>::iterator
MathVector<kNumberOfElements, ScalarType>::end(void) noexcept {
  return elements_m.end();
}
template <UnsignedIndex_t kNumberOfElements, class ScalarType>
typename MathVector<kNumberOfElements, ScalarType>::const_iterator
MathVector<kNumberOfElements, ScalarType>::end(void) const noexcept {
  return this->cend();
}
template <UnsignedIndex_t kNumberOfElements, class ScalarType>
typename MathVector<kNumberOfElements, ScalarType>::const_iterator
MathVector<kNumberOfElements, ScalarType>::cend(void) const noexcept {
  return elements_m.cend();
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
LargeOffsetIndex_t MathVector<kNumberOfElements, ScalarType>::getSerializedSize(
    void) const {
  return static_cast<LargeOffsetIndex_t>(this->size() * sizeof(ScalarType));
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
void MathVector<kNumberOfElements, ScalarType>::serialize(
    ByteBuffer* a_buffer) const {
  a_buffer->pack(elements_m.data(), this->size());
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
void MathVector<kNumberOfElements, ScalarType>::unpackSerialized(
    ByteBuffer* a_buffer) {
  a_buffer->unpack(elements_m.data(), this->size());
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType>::MathVector(const ScalarType* a_vec) {
  for (UnsignedIndex_t ind = 0; ind < kNumberOfElements; ++ind) {
    elements_m[ind] = a_vec[ind];
  }
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType>::MathVector(
    const ScalarType a_constant) {
  std::fill(this->begin(), this->end(), a_constant);
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
ScalarType operator*(const MathVector<kNumberOfElements, ScalarType>& a_vec_0,
                     const MathVector<kNumberOfElements, ScalarType>& a_vec_1) {
  ScalarType dot_product = static_cast<ScalarType>(0.0);
  for (UnsignedIndex_t ind = 0; ind < a_vec_0.size(); ++ind) {
    dot_product += a_vec_0[ind] * a_vec_1[ind];
  }
  return dot_product;
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType> operator+(
    const MathVector<kNumberOfElements, ScalarType>& a_vec_0,
    const MathVector<kNumberOfElements, ScalarType>& a_vec_1) {
  MathVector<kNumberOfElements, ScalarType> sum;
  for (UnsignedIndex_t n = 0; n < a_vec_0.size(); ++n) {
    sum[n] = a_vec_0[n] + a_vec_1[n];
  }
  return sum;
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
const MathVector<kNumberOfElements, ScalarType> operator-(
    const MathVector<kNumberOfElements, ScalarType>& a_vec_0,
    const MathVector<kNumberOfElements, ScalarType>& a_vec_1) {
  MathVector<kNumberOfElements, ScalarType> difference;
  for (UnsignedIndex_t n = 0; n < a_vec_0.size(); ++n) {
    difference[n] = a_vec_0[n] - a_vec_1[n];
  }
  return difference;
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType> operator*(
    const ScalarType a_scalar,
    const MathVector<kNumberOfElements, ScalarType>& a_vec) {
  auto vec_to_return = a_vec;
  vec_to_return *= a_scalar;
  return vec_to_return;
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType> operator*(
    const MathVector<kNumberOfElements, ScalarType>& a_vec,
    const ScalarType a_scalar) {
  return a_scalar * a_vec;
}

template <UnsignedIndex_t kNumberOfElements, class ScalarType>
MathVector<kNumberOfElements, ScalarType> operator/(
    const MathVector<kNumberOfElements, ScalarType>& a_vec,
    const ScalarType a_scalar) {
  MathVector<kNumberOfElements, ScalarType> vec_to_return = a_vec;
  vec_to_return /= a_scalar;
  return vec_to_return;
}
}  // namespace IRL

#endif  // IRL_GEOMETRY_GENERAL_MATH_VECTOR_TPP_
