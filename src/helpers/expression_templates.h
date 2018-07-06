// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_HELPERS_EXPRESSION_TEMPLATES_H_
#define SRC_HELPERS_EXPRESSION_TEMPLATES_H_

#include <float.h>
#include <type_traits>
#include <utility>

#include <cassert>
#include "src/parameters/defined_types.h"

namespace IRL {

template <class E>
struct Expr {
  operator const E&() const;
};

// Class from Discovering Modern C++, Peter Gottschling, 2016
template <typename V1, typename V2>
class vector_sum : public Expr<vector_sum<V1, V2>> {
  using self = vector_sum;

 public:
  using value_type = typename std::common_type<typename V1::value_type,
                                               typename V2::value_type>::type;

  friend UnsignedIndex_t size(const self& x) { return size(x.v1_m); }

  vector_sum(const V1& a_v1, const V2& a_v2);

  value_type operator[](const UnsignedIndex_t i) const;

 private:
  void check_index(const UnsignedIndex_t i) const;
  const V1& v1_m;
  const V2& v2_m;
};

template <typename V1, typename V2>
enable_if_t<std::is_base_of<Expr<V1>, V1>::value &&
                std::is_base_of<Expr<V2>, V2>::value,
            vector_sum<V1, V2>> inline
operator+(const Expr<V1>& x, const Expr<V2>& y);

template <typename V1, typename V2>
class vector_subtract : public Expr<vector_subtract<V1, V2>> {
  using self = vector_subtract;

 public:
  using value_type = typename std::common_type<typename V1::value_type,
                                               typename V2::value_type>::type;

  friend UnsignedIndex_t size(const self& x) { return size(x.v1_m); }

  vector_subtract(const V1& a_v1, const V2& a_v2);

  value_type operator[](const UnsignedIndex_t i) const;

 private:
  void check_index(const UnsignedIndex_t i) const;
  const V1& v1_m;
  const V2& v2_m;
};

template <typename V1, typename V2>
enable_if_t<std::is_base_of<Expr<V1>, V1>::value &&
                std::is_base_of<Expr<V2>, V2>::value,
            vector_subtract<V1, V2>> inline
operator-(const Expr<V1>& x, const Expr<V2>& y);

template <typename V1>
class vector_scale : public Expr<vector_scale<V1>> {
  using self = vector_scale;

 public:
  using value_type = typename V1::value_type;

  friend UnsignedIndex_t size(const self& x) { return size(x.v1_m); }

  vector_scale(const double a_scalar, const V1& a_v1);

  value_type operator[](const UnsignedIndex_t i) const;

 private:
  void check_index(const UnsignedIndex_t i) const;
  const double scalar_m;
  const V1& v1_m;
};

template <typename V1>
enable_if_t<std::is_base_of<Expr<V1>, V1>::value,
            vector_scale<V1>> inline
operator*(const double a_scalar, const Expr<V1>& x);

template <typename V1>
enable_if_t<std::is_base_of<Expr<V1>, V1>::value,
            vector_scale<V1>> inline
operator*(const Expr<V1>& x, const double a_scalar);


template <typename V1, typename V2>
class vector_cross_product : public Expr<vector_cross_product<V1, V2>> {
  using self = vector_cross_product;

 public:
  using value_type = typename std::common_type<typename V1::value_type,
                                               typename V2::value_type>::type;

  friend UnsignedIndex_t size(const self& x) { return size(x.v1_m); }

  vector_cross_product(const V1& a_v1, const V2& a_v2);

  value_type operator[](const UnsignedIndex_t i) const;

 private:
  void check_index(const UnsignedIndex_t i) const;
  const V1& v1_m;
  const V2& v2_m;
};

template <typename V1, typename V2>
enable_if_t<std::is_base_of<Expr<V1>, V1>::value &&
                std::is_base_of<Expr<V2>, V2>::value,
				vector_cross_product<V1, V2>> inline
crossProduct(const Expr<V1>& x, const Expr<V2>& y);

}  // namespace IRL

#include "src/helpers/expression_templates.tpp"

#endif  // SRC_HELPERS_EXPRESSION_TEMPLATES_H_
