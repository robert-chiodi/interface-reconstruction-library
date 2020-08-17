// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_MATH_VECTOR_TPP_
#define SRC_GEOMETRY_GENERAL_MATH_VECTOR_TPP_

#include <algorithm>

namespace IRL {

template <UnsignedIndex_t kNumberOfElements>
constexpr MathVector<kNumberOfElements>::MathVector(const double a_element_0,
                                                    const double a_element_1,
                                                    const double a_element_2)
    : elements_m{a_element_0, a_element_1, a_element_2} {}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements>
MathVector<kNumberOfElements>::fromRawDoublePointer(const double* a_vec) {
  return MathVector(a_vec);
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements> MathVector<kNumberOfElements>::fromScalarConstant(
    const double a_constant) {
  return MathVector(a_constant);
}

template <UnsignedIndex_t kNumberOfElements>
template <class E>
inline MathVector<kNumberOfElements>::MathVector(const Expr<E>& a_expr) {
  const E& expr(a_expr);
  for(UnsignedIndex_t n = 0; n < kNumberOfElements; ++n){
    elements_m[n] = expr[n];
  }
}

template <UnsignedIndex_t kNumberOfElements>
template <class E>
inline MathVector<kNumberOfElements>::MathVector(Expr<E>&& a_expr) {
  const E& expr(a_expr);
  for(UnsignedIndex_t n = 0; n < kNumberOfElements; ++n){
    elements_m[n] = expr[n];
  }
}

template <UnsignedIndex_t kNumberOfElements>
template <class E>
inline MathVector<kNumberOfElements>& MathVector<kNumberOfElements>::operator=(const Expr<E>& a_expr) {
  const E& expr(a_expr);
  for(UnsignedIndex_t n = 0; n < kNumberOfElements; ++n){
    elements_m[n] = expr[n];
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements>
template <class E>
inline MathVector<kNumberOfElements>& MathVector<kNumberOfElements>::operator=(Expr<E>&& a_expr) {
  const E& expr(a_expr);
  for(UnsignedIndex_t n = 0; n < kNumberOfElements; ++n){
    elements_m[n] = expr[n];
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements>
constexpr UnsignedIndex_t MathVector<kNumberOfElements>::size(void) {
  return kNumberOfElements;
}

template <UnsignedIndex_t kNumberOfElements>
double& MathVector<kNumberOfElements>::operator[](const UnsignedIndex_t a_d) {
  assert(a_d < this->size());
  return elements_m[a_d];
}

template <UnsignedIndex_t kNumberOfElements>
const double& MathVector<kNumberOfElements>::operator[](
    const UnsignedIndex_t a_d) const {
  assert(a_d < this->size());
  return elements_m[a_d];
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements> MathVector<kNumberOfElements>::operator-(
    void) const {
  MathVector vector_to_return = (*this);
  for (auto& element : (*this)) {
    element = -element;
  }
  return vector_to_return;
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements>& MathVector<kNumberOfElements>::operator+=(
    const MathVector& a_other_vector) {
  assert(a_other_vector.size() == this->size());
  for (UnsignedIndex_t ind = 0; ind < a_other_vector.size(); ++ind) {
    (*this)[ind] += a_other_vector[ind];
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements>& MathVector<kNumberOfElements>::operator-=(
    const MathVector& a_other_vector) {
  assert(a_other_vector.size() == this->size());
  for (UnsignedIndex_t ind = 0; ind < a_other_vector.size(); ++ind) {
    (*this)[ind] -= a_other_vector[ind];
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements>& MathVector<kNumberOfElements>::operator+=(
    const double a_scalar) {
  for (auto& element : (*this)) {
    element += a_scalar;
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements>& MathVector<kNumberOfElements>::operator-=(
    const double a_scalar) {
  for (auto& element : (*this)) {
    element -= a_scalar;
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements>& MathVector<kNumberOfElements>::operator*=(
    const double a_scalar) {
  for (auto& element : (*this)) {
    element *= a_scalar;
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements>& MathVector<kNumberOfElements>::operator/=(
    const double a_scalar) {
  for (auto& element : (*this)) {
    element /= a_scalar;
  }
  return (*this);
}

template <UnsignedIndex_t kNumberOfElements>
UnsignedIndex_t MathVector<kNumberOfElements>::calculateIndexOfLargestMagnitude(
    void) const {
  UnsignedIndex_t index_of_max_magnitude;
  double max_value = -DBL_MAX;
  double abs_value_of_element;
  for (UnsignedIndex_t ind = 0; ind < this->size(); ++ind) {
    abs_value_of_element = std::fabs(elements_m[ind]);
    if (abs_value_of_element > max_value) {
      max_value = abs_value_of_element;
      index_of_max_magnitude = ind;
    }
  }
  return index_of_max_magnitude;
}

template <UnsignedIndex_t kNumberOfElements>
double MathVector<kNumberOfElements>::calculateMagnitude(void) const {
  return magnitude((*this));
}

template <UnsignedIndex_t kNumberOfElements>
typename MathVector<kNumberOfElements>::iterator
MathVector<kNumberOfElements>::begin(void) noexcept {
  return elements_m.begin();
}
template <UnsignedIndex_t kNumberOfElements>
typename MathVector<kNumberOfElements>::const_iterator
MathVector<kNumberOfElements>::begin(void) const noexcept {
  return this->cbegin();
}
template <UnsignedIndex_t kNumberOfElements>
typename MathVector<kNumberOfElements>::const_iterator
MathVector<kNumberOfElements>::cbegin(void) const noexcept {
  return elements_m.cbegin();
}
template <UnsignedIndex_t kNumberOfElements>
typename MathVector<kNumberOfElements>::iterator
MathVector<kNumberOfElements>::end(void) noexcept {
  return elements_m.end();
}
template <UnsignedIndex_t kNumberOfElements>
typename MathVector<kNumberOfElements>::const_iterator
MathVector<kNumberOfElements>::end(void) const noexcept {
  return this->cend();
}
template <UnsignedIndex_t kNumberOfElements>
typename MathVector<kNumberOfElements>::const_iterator
MathVector<kNumberOfElements>::cend(void) const noexcept {
  return elements_m.cend();
}

template <UnsignedIndex_t kNumberOfElements>
LargeOffsetIndex_t MathVector<kNumberOfElements>::getSerializedSize(
    void) const {
  return static_cast<LargeOffsetIndex_t>(this->size() * sizeof(double));
}

template <UnsignedIndex_t kNumberOfElements>
void MathVector<kNumberOfElements>::serialize(ByteBuffer* a_buffer) const {
  a_buffer->pack(elements_m.data(), this->size());
}

template <UnsignedIndex_t kNumberOfElements>
void MathVector<kNumberOfElements>::unpackSerialized(ByteBuffer* a_buffer) {
  a_buffer->unpack(elements_m.data(), this->size());
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements>::MathVector(const double* a_vec) {
  for (UnsignedIndex_t ind = 0; ind < kNumberOfElements; ++ind) {
    elements_m[ind] = a_vec[ind];
  }
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements>::MathVector(const double a_constant) {
  std::fill(this->begin(), this->end(), a_constant);
}

template <UnsignedIndex_t kNumberOfElements>
double operator*(const MathVector<kNumberOfElements>& a_vec_0,
                 const MathVector<kNumberOfElements>& a_vec_1) {
  double dot_product = 0.0;
  for (UnsignedIndex_t ind = 0; ind < a_vec_0.size(); ++ind) {
    dot_product += a_vec_0[ind] * a_vec_1[ind];
  }
  return dot_product;
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements> operator+(
    const MathVector<kNumberOfElements>& a_vec_0,
    const MathVector<kNumberOfElements>& a_vec_1) {
  MathVector<kNumberOfElements> sum;
  for (UnsignedIndex_t n = 0; n < a_vec_0.size(); ++n) {
    sum[n] = a_vec_0[n] + a_vec_1[n];
  }
  return sum;
}

template <UnsignedIndex_t kNumberOfElements>
const MathVector<kNumberOfElements> operator-(
    const MathVector<kNumberOfElements>& a_vec_0,
    const MathVector<kNumberOfElements>& a_vec_1) {
  MathVector<kNumberOfElements> difference;
  for (UnsignedIndex_t n = 0; n < a_vec_0.size(); ++n) {
    difference[n] = a_vec_0[n] - a_vec_1[n];
  }
  return difference;
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements> operator*(
    const double a_double, const MathVector<kNumberOfElements>& a_vec) {
  auto vec_to_return = a_vec;
  vec_to_return *= a_double;
  return vec_to_return;
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements> operator*(
    const MathVector<kNumberOfElements>& a_vec, const double a_double) {
  return a_double * a_vec;
}

template <UnsignedIndex_t kNumberOfElements>
MathVector<kNumberOfElements> operator/(
    const MathVector<kNumberOfElements>& a_vec, const double a_double) {
  auto vec_to_return = a_vec;
  vec_to_return /= a_double;
  return vec_to_return;
}
}  // namespace IRL

#endif  // SRC_GEOMETRY_GENERAL_MATH_VECTOR_TPP_
