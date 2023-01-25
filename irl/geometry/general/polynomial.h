// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi  <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_POLYNOMIAL_H_
#define IRL_GEOMETRY_GENERAL_POLYNOMIAL_H_

#include <algorithm>

namespace IRL {

template <UnsignedIndex_t kNCoefficients, class ScalarType>
class PolynomialBase {
 public:
  PolynomialBase(void) {
    std::fill(coefficients_m.begin(), coefficients_m.end(),
              static_cast<ScalarType>(0));
  }

  PolynomialBase(const std::array<ScalarType, kNCoefficients>& a_coefficients) {
    coefficients_m = a_coefficients;
  }

  ScalarType& operator[](const UnsignedIndex_t a_index) {
    assert(a_index < coefficients_m.size());
    return coefficients_m[a_index];
  }
  ScalarType operator[](const UnsignedIndex_t a_index) const {
    assert(a_index < coefficients_m.size());
    return coefficients_m[a_index];
  }
  auto begin(void) { return coefficients_m.begin(); }
  auto begin(void) const { return this->cbegin(); }
  auto cbegin(void) const { return coefficients_m.cbegin(); }
  auto end(void) { return coefficients_m.end(); }
  auto end(void) const { return this->cend(); }
  auto cend(void) const { return coefficients_m.cend(); }

  ~PolynomialBase(void) = default;

 private:
  std::array<ScalarType, kNCoefficients> coefficients_m;
};

template <UnsignedIndex_t kNCoefficients>
using Polynomial = PolynomialBase<kNCoefficients, double>;

}  // namespace IRL

#include "irl/geometry/general/polynomial.tpp"

#endif  // IRL_GEOMETRY_GENERAL_POLYNOMIAL_H_
