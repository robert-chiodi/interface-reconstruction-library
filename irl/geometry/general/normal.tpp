// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_NORMAL_TPP_
#define IRL_GEOMETRY_GENERAL_NORMAL_TPP_

namespace IRL {

template <class ScalarType>
inline constexpr NormalBase<ScalarType>::NormalBase(void)
    : normal_m{static_cast<ScalarType>(0), static_cast<ScalarType>(0),
               static_cast<ScalarType>(0)} {}

template <class ScalarType>
inline constexpr NormalBase<ScalarType>::NormalBase(
    const UnsignedIndex_t a_normal_x, const UnsignedIndex_t a_normal_y,
    const UnsignedIndex_t a_normal_z)
    : normal_m{static_cast<ScalarType>(a_normal_x),
               static_cast<ScalarType>(a_normal_y),
               static_cast<ScalarType>(a_normal_z)} {}

template <class ScalarType>
inline constexpr NormalBase<ScalarType>::NormalBase(const int a_normal_x,
                                                    const int a_normal_y,
                                                    const int a_normal_z)
    : normal_m{static_cast<ScalarType>(a_normal_x),
               static_cast<ScalarType>(a_normal_y),
               static_cast<ScalarType>(a_normal_z)} {}

template <class ScalarType>
inline constexpr NormalBase<ScalarType>::NormalBase(const double a_normal_x,
                                                    const double a_normal_y,
                                                    const double a_normal_z)
    : normal_m{a_normal_x, a_normal_y, a_normal_z} {
  static_assert(
      std::is_same<ScalarType, double>::value,
      "Trying to construct Normal in quad precision with doubles. This "
      "is forbidden to avoid loosing quad precision.");
}

template <class ScalarType>
inline constexpr NormalBase<ScalarType>::NormalBase(const Quad_t a_normal_x,
                                                    const Quad_t a_normal_y,
                                                    const Quad_t a_normal_z)
    : normal_m{static_cast<ScalarType>(a_normal_x),
               static_cast<ScalarType>(a_normal_y),
               static_cast<ScalarType>(a_normal_z)} {}

template <class ScalarType>
inline NormalBase<ScalarType> NormalBase<ScalarType>::normalized(
    const ScalarType a_normal_x, const ScalarType a_normal_y,
    const ScalarType a_normal_z) {
  auto normal_to_normalize =
      NormalBase<ScalarType>(a_normal_x, a_normal_y, a_normal_z);
  normal_to_normalize.normalize();
  return normal_to_normalize;
}

template <class ScalarType>
inline NormalBase<ScalarType> NormalBase<ScalarType>::fromRawDoublePointer(
    const ScalarType* a_normal) {
  return NormalBase<ScalarType>(a_normal);
}

template <class ScalarType>
inline NormalBase<ScalarType>
NormalBase<ScalarType>::fromRawDoublePointerNormalized(
    const ScalarType* a_normal) {
  auto normal_to_normalize = NormalBase<ScalarType>(a_normal);
  normal_to_normalize.normalize();
  return normal_to_normalize;
}

template <class ScalarType>
inline NormalBase<ScalarType> NormalBase<ScalarType>::fromPt(
    const PtBase<ScalarType>& a_pt) {
  return NormalBase<ScalarType>(a_pt);
}

template <class ScalarType>
inline NormalBase<ScalarType> NormalBase<ScalarType>::fromPtNormalized(
    const PtBase<ScalarType>& a_pt) {
  auto normal_to_normalize = NormalBase<ScalarType>(a_pt);
  normal_to_normalize.normalize();
  return normal_to_normalize;
}

template <class ScalarType>
inline PtBase<ScalarType> NormalBase<ScalarType>::toPt(
    const NormalBase<ScalarType>& a_normal) {
  return PtBase<ScalarType>(a_normal[0], a_normal[1], a_normal[2]);
}

template <class ScalarType>
inline const NormalBase<double> NormalBase<ScalarType>::toDoubleNormal(
    void) const {
  if constexpr (std::is_same<ScalarType, double>::value) {
    return (*this);
  } else {
    return NormalBase<double>(static_cast<double>(normal_m[0]),
                              static_cast<double>(normal_m[1]),
                              static_cast<double>(normal_m[2]));
  }
}

template <class ScalarType>
inline const NormalBase<Quad_t> NormalBase<ScalarType>::toQuadNormal(
    void) const {
  if constexpr (std::is_same<ScalarType, Quad_t>::value) {
    return (*this);
  } else {
    return NormalBase<Quad_t>(static_cast<Quad_t>(normal_m[0]),
                              static_cast<Quad_t>(normal_m[1]),
                              static_cast<Quad_t>(normal_m[2]));
  }
}

template <class ScalarType>
inline NormalBase<ScalarType> NormalBase<ScalarType>::fromScalarConstant(
    const ScalarType a_constant) {
  return NormalBase<ScalarType>(a_constant);
}

template <class ScalarType>
template <class E>
inline NormalBase<ScalarType>::NormalBase(const Expr<E>& a_expr)
    : normal_m(a_expr) {}

template <class ScalarType>
template <class E>
inline NormalBase<ScalarType>::NormalBase(Expr<E>&& a_expr)
    : normal_m(std::move(a_expr)) {}

template <class ScalarType>
template <class E>
inline NormalBase<ScalarType>& NormalBase<ScalarType>::operator=(
    const Expr<E>& a_expr) {
  normal_m = a_expr;
  return (*this);
}

template <class ScalarType>
template <class E>
inline NormalBase<ScalarType>& NormalBase<ScalarType>::operator=(
    Expr<E>&& a_expr) {
  normal_m = std::move(a_expr);
  return (*this);
}

template <class ScalarType>
inline ScalarType& NormalBase<ScalarType>::operator[](
    const UnsignedIndex_t a_d) {
  assert(a_d < 3);
  return normal_m[a_d];
}

template <class ScalarType>
inline const ScalarType& NormalBase<ScalarType>::operator[](
    const UnsignedIndex_t a_d) const {
  assert(a_d < 3);
  return normal_m[a_d];
}

template <class ScalarType>
inline constexpr UnsignedIndex_t NormalBase<ScalarType>::size(void) {
  return 3;
}

template <class ScalarType>
inline NormalBase<ScalarType> NormalBase<ScalarType>::operator-(void) const {
  return {-normal_m[0], -normal_m[1], -normal_m[2]};
}

template <class ScalarType>
inline NormalBase<ScalarType>& NormalBase<ScalarType>::operator+=(
    const NormalBase<ScalarType>& a_rhs) {
  normal_m[0] += a_rhs[0];
  normal_m[1] += a_rhs[1];
  normal_m[2] += a_rhs[2];
  return (*this);
}

template <class ScalarType>
inline NormalBase<ScalarType>& NormalBase<ScalarType>::operator/=(
    const double a_value) {
  static_assert(std::is_same<ScalarType, double>::value,
                "Trying to divide Normal in quad precision with double. This "
                "is forbidden to avoid loosing quad precision.");
  assert(a_value != 0.0);
  for (auto& elem : normal_m) {
    elem /= a_value;
  }
  return (*this);
}

template <class ScalarType>
inline NormalBase<ScalarType>& NormalBase<ScalarType>::operator/=(
    const Quad_t a_value) {
  assert(a_value != 0.0q);
  for (auto& elem : normal_m) {
    elem /= static_cast<ScalarType>(a_value);
  }
  return (*this);
}

template <class ScalarType>
inline NormalBase<ScalarType>& NormalBase<ScalarType>::operator*=(
    const double a_value) {
  static_assert(std::is_same<ScalarType, double>::value,
                "Trying to multiply Normal in quad precision with double. This "
                "is forbidden to avoid loosing quad precision.");
  for (auto& elem : normal_m) {
    elem *= a_value;
  }
  return (*this);
}

template <class ScalarType>
inline NormalBase<ScalarType>& NormalBase<ScalarType>::operator*=(
    const Quad_t a_value) {
  for (auto& elem : normal_m) {
    elem *= static_cast<ScalarType>(a_value);
  }
  return (*this);
}

template <class ScalarType>
inline NormalBase<ScalarType>& NormalBase<ScalarType>::operator=(
    const double a_value) {
  static_assert(std::is_same<ScalarType, double>::value,
                "Trying to equate Normal in quad precision with double. This "
                "is forbidden to avoid loosing quad precision.");
  for (auto& elem : normal_m) {
    elem = a_value;
  }
  return (*this);
}

template <class ScalarType>
inline NormalBase<ScalarType>& NormalBase<ScalarType>::operator=(
    const Quad_t a_value) {
  for (auto& elem : normal_m) {
    elem = static_cast<ScalarType>(a_value);
  }
  return (*this);
}

template <class ScalarType>
inline bool NormalBase<ScalarType>::operator==(
    const NormalBase<ScalarType>& a_normal) const {
  return (*this) * a_normal >
         static_cast<ScalarType>(global_constants::SAME_VEC);
}

template <class ScalarType>
inline bool NormalBase<ScalarType>::operator!=(
    const NormalBase<ScalarType>& a_normal) const {
  return (*this) * a_normal <
         static_cast<ScalarType>(global_constants::SAME_VEC);
}

template <class ScalarType>
inline UnsignedIndex_t NormalBase<ScalarType>::maxMagnitudeComponent(
    void) const {
  UnsignedIndex_t max = fabs((*this)[0]) > fabs((*this)[1]) ? 0 : 1;
  max = fabs((*this)[max]) > fabs((*this)[2]) ? max : 2;
  return max;
}

template <class ScalarType>
inline ScalarType NormalBase<ScalarType>::calculateMagnitude(void) const {
  return magnitude((*this));
}

template <class ScalarType>
inline ScalarType NormalBase<ScalarType>::calculateSquaredMagnitude(
    void) const {
  return squaredMagnitude((*this));
}

template <>
inline void NormalBase<double>::normalize(void) {
  // const double inv_magnitude = 1.0 / safelyTiny(this->calculateMagnitude());
  const double inv_magnitude =
      invsqrt(safelyTiny(this->calculateSquaredMagnitude()));
  for (auto& elem : normal_m) {
    elem *= inv_magnitude;
  }
}

template <>
inline void NormalBase<Quad_t>::normalize(void) {
  const Quad_t inv_magnitude = 1.0q / safelyTiny(this->calculateMagnitude());
  for (auto& elem : normal_m) {
    elem *= inv_magnitude;
  }
}

template <class ScalarType>
inline typename NormalBase<ScalarType>::iterator NormalBase<ScalarType>::begin(
    void) noexcept {
  return normal_m.begin();
}
template <class ScalarType>
inline typename NormalBase<ScalarType>::const_iterator
NormalBase<ScalarType>::begin(void) const noexcept {
  return this->cbegin();
}
template <class ScalarType>
inline typename NormalBase<ScalarType>::const_iterator
NormalBase<ScalarType>::cbegin(void) const noexcept {
  return normal_m.cbegin();
}
template <class ScalarType>
inline typename NormalBase<ScalarType>::iterator NormalBase<ScalarType>::end(
    void) noexcept {
  return normal_m.end();
}
template <class ScalarType>
inline typename NormalBase<ScalarType>::const_iterator
NormalBase<ScalarType>::end(void) const noexcept {
  return this->cend();
}
template <class ScalarType>
inline typename NormalBase<ScalarType>::const_iterator
NormalBase<ScalarType>::cend(void) const noexcept {
  return normal_m.cend();
}

template <class ScalarType>
inline LargeOffsetIndex_t NormalBase<ScalarType>::getSerializedSize(
    void) const {
  return normal_m.getSerializedSize();
}

template <class ScalarType>
inline void NormalBase<ScalarType>::serialize(ByteBuffer* a_buffer) const {
  normal_m.serialize(a_buffer);
}

template <class ScalarType>
inline void NormalBase<ScalarType>::unpackSerialized(ByteBuffer* a_buffer) {
  normal_m.unpackSerialized(a_buffer);
}

template <class ScalarType>
inline constexpr NormalBase<ScalarType>::NormalBase(const double* a_normal)
    : normal_m{a_normal[0], a_normal[1], a_normal[2]} {
  static_assert(
      std::is_same<ScalarType, double>::value,
      "Trying to construct Normal in quad precision with doubles. This "
      "is forbidden to avoid loosing quad precision.");
}

template <class ScalarType>
inline constexpr NormalBase<ScalarType>::NormalBase(const Quad_t* a_normal)
    : normal_m{static_cast<ScalarType>(a_normal[0]),
               static_cast<ScalarType>(a_normal[1]),
               static_cast<ScalarType>(a_normal[2])} {}

template <class ScalarType>
inline NormalBase<ScalarType>::NormalBase(const PtBase<ScalarType>& a_pt)
    : normal_m{a_pt[0], a_pt[1], a_pt[2]} {}

template <class ScalarType>
inline NormalBase<ScalarType>::NormalBase(const double a_constant)
    : normal_m{a_constant, a_constant, a_constant} {
  static_assert(
      std::is_same<ScalarType, double>::value,
      "Trying to construct Normal in quad precision with double. This "
      "is forbidden to avoid loosing quad precision.");
}

template <class ScalarType>
inline NormalBase<ScalarType>::NormalBase(const Quad_t a_constant)
    : normal_m{static_cast<ScalarType>(a_constant),
               static_cast<ScalarType>(a_constant),
               static_cast<ScalarType>(a_constant)} {}

template <class ScalarType>
__attribute__((const)) inline ScalarType operator*(
    const NormalBase<ScalarType>& a_normal_0,
    const NormalBase<ScalarType>& a_normal_1) {
  return dotProduct(a_normal_0, a_normal_1);
}

template <class ScalarType>
__attribute__((const)) inline ScalarType operator*(
    const PtBase<ScalarType>& a_pt, const NormalBase<ScalarType>& a_normal) {
  return dotProduct(a_pt, a_normal);
}

template <class ScalarType>
__attribute__((const)) inline ScalarType operator*(
    const NormalBase<ScalarType>& a_normal, const PtBase<ScalarType>& a_pt) {
  return a_pt * a_normal;
}

template <class ScalarType>
__attribute__((const)) inline NormalBase<ScalarType> operator*(
    const ScalarType a_double, const NormalBase<ScalarType>& a_normal) {
  return NormalBase<ScalarType>(a_normal[0] * a_double, a_normal[1] * a_double,
                                a_normal[2] * a_double);
}

template <class ScalarType>
__attribute__((const)) inline NormalBase<ScalarType> operator*(
    const NormalBase<ScalarType>& a_normal, const ScalarType a_double) {
  return a_double * a_normal;
}

template <class ScalarType>
__attribute__((const)) inline NormalBase<ScalarType> operator/(
    const NormalBase<ScalarType>& a_normal, const ScalarType a_double) {
  return NormalBase<ScalarType>(a_normal[0] / a_double, a_normal[1] / a_double,
                                a_normal[2] / a_double);
}

template <class ScalarType>
inline std::ostream& operator<<(std::ostream& out,
                                const NormalBase<ScalarType>& a_normal) {
  out << std::setprecision(15);
  out << "( " << a_normal[0];
  out << ", " << a_normal[1];
  out << ", " << a_normal[2];
  out << " )";
  return out;
}

}  // namespace IRL

#endif  // IRL_GEOMETRY_GENERAL_NORMAL_TPP_
