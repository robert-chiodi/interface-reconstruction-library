// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi  <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_POLYNOMIAL_TPP_
#define IRL_GEOMETRY_GENERAL_POLYNOMIAL_TPP_

namespace IRL {

template <UnsignedIndex_t kNCoefficients, class ScalarType>
inline void PolynomialBase<kNCoefficients, ScalarType>::serialize(
    ByteBuffer* a_buffer) const {
  a_buffer->pack(coefficients_m.data(), kNCoefficients);
}

template <UnsignedIndex_t kNCoefficients, class ScalarType>
inline void PolynomialBase<kNCoefficients, ScalarType>::unpackSerialized(
    ByteBuffer* a_buffer) {
  a_buffer->unpack(coefficients_m.data(), kNCoefficients);
}
}  // namespace IRL

#endif  // IRL_GEOMETRY_GENERAL_POLYNOMIAL_TPP_
