// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_HELPERS_EXPRESSION_TEMPLATES_TPP_
#define SRC_HELPERS_EXPRESSION_TEMPLATES_TPP_

namespace IRL {

template <class E>
Expr<E>::operator const E&() const {
  return static_cast<const E&>(*this);
}

template <typename V1, typename V2>
vector_sum<V1, V2>::vector_sum(const V1& a_v1, const V2& a_v2)
    : v1_m(a_v1), v2_m(a_v2) {
  assert(size(v1_m) == size(v2_m));
}

template <typename V1, typename V2>
typename vector_sum<V1, V2>::value_type vector_sum<V1, V2>::operator[](
    const UnsignedIndex_t i) const {
  check_index(i);
  return v1_m[i] + v2_m[i];
}

template <typename V1, typename V2>
void vector_sum<V1, V2>::check_index(const UnsignedIndex_t i) const {
  assert(i < size(v1_m));
}

template <typename V1, typename V2>
enable_if_t<std::is_base_of<Expr<V1>, V1>::value &&
                std::is_base_of<Expr<V2>, V2>::value,
            vector_sum<V1, V2>> inline
operator+(const Expr<V1>& x, const Expr<V2>& y) {
  return vector_sum<V1, V2>(x, y);
}

template <typename V1, typename V2>
vector_subtract<V1, V2>::vector_subtract(const V1& a_v1, const V2& a_v2)
    : v1_m(a_v1), v2_m(a_v2) {
  assert(size(v1_m) == size(v2_m));
}

template <typename V1, typename V2>
typename vector_subtract<V1, V2>::value_type vector_subtract<V1, V2>::operator[](
    const UnsignedIndex_t i) const {
  check_index(i);
  return v1_m[i] - v2_m[i];
}

template <typename V1, typename V2>
void vector_subtract<V1, V2>::check_index(const UnsignedIndex_t i) const {
  assert(i < size(v1_m));
}

template <typename V1, typename V2>
enable_if_t<std::is_base_of<Expr<V1>, V1>::value &&
                std::is_base_of<Expr<V2>, V2>::value,
            vector_subtract<V1, V2>> inline
operator-(const Expr<V1>& x, const Expr<V2>& y) {
  return vector_subtract<V1, V2>(x, y);
}


template <typename V1>
vector_scale<V1>::vector_scale(const double a_scalar, const V1& a_v1)
    : scalar_m(a_scalar), v1_m(a_v1) {}

template <typename V1>
typename vector_scale<V1>::value_type vector_scale<V1>::operator[](
    const UnsignedIndex_t i) const {
  check_index(i);
  return scalar_m*v1_m[i];
}

template <typename V1>
void vector_scale<V1>::check_index(const UnsignedIndex_t i) const {
  assert(i < size(v1_m));
}

template <typename V1>
enable_if_t<std::is_base_of<Expr<V1>, V1>::value,
            vector_scale<V1>> inline
operator*(const double a_scalar, const Expr<V1>& x) {
  return vector_scale<V1>(a_scalar, x);
}

template <typename V1>
enable_if_t<std::is_base_of<Expr<V1>, V1>::value,
            vector_scale<V1>> inline
operator*(const Expr<V1>& x, const double a_scalar) {
  return vector_scale<V1>(a_scalar, x);
}

template <typename V1, typename V2>
vector_cross_product<V1, V2>::vector_cross_product(const V1& a_v1, const V2& a_v2)
    : v1_m(a_v1), v2_m(a_v2) {
  assert(size(v1_m) == size(v2_m));
}

template <typename V1, typename V2>
typename vector_cross_product<V1, V2>::value_type vector_cross_product<V1, V2>::operator[](
    const UnsignedIndex_t i) const {
  check_index(i);
  assert(size(v1_m) == 3);
  assert(size(v2_m) == 3);
  const auto j = (i+1)%3;
  const auto k = (i+2)%3;
  return v1_m[j] * v2_m[k] - v1_m[k] * v2_m[j];

//  1: v1_m[1] * v2_m[2] - v1_m[2] * v2_m[1];
//  2: v1_m[2] * v2_m[0] - v1_m[0] * v2_m[2];
//  3: v1_m[0] * v2_m[1] - v1_m[1] * v2_m[0];
}

template <typename V1, typename V2>
void vector_cross_product<V1, V2>::check_index(const UnsignedIndex_t i) const {
  assert(i < size(v1_m));
}

template <typename V1, typename V2>
enable_if_t<std::is_base_of<Expr<V1>, V1>::value &&
                std::is_base_of<Expr<V2>, V2>::value,
				vector_cross_product<V1, V2>> inline
  crossProduct(const Expr<V1>& x, const Expr<V2>& y) {
  return vector_cross_product<V1, V2>(x, y);
}
	
}  // namespace IRL

#endif  // SRC_HELPERS_EXPRESSION_TEMPLATES_TPP_
