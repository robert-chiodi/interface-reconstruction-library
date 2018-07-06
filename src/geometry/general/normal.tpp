// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_NORMAL_TPP_
#define SRC_GEOMETRY_GENERAL_NORMAL_TPP_

namespace IRL {

inline constexpr Normal::Normal(void) : normal_m{0.0, 0.0, 0.0} {}

inline constexpr Normal::Normal(const double a_normal_x,
                                const double a_normal_y,
                                const double a_normal_z)
    : normal_m{a_normal_x, a_normal_y, a_normal_z} {}

inline Normal Normal::normalized(const double a_normal_x,
                                 const double a_normal_y,
                                 const double a_normal_z) {
  auto normal_to_normalize = Normal(a_normal_x, a_normal_y, a_normal_z);
  normal_to_normalize.normalize();
  return normal_to_normalize;
}

inline Normal Normal::fromRawDoublePointer(const double* a_normal) {
  return Normal(a_normal);
}

inline Normal Normal::fromRawDoublePointerNormalized(const double* a_normal) {
  auto normal_to_normalize = Normal(a_normal);
  normal_to_normalize.normalize();
  return normal_to_normalize;
}

inline Normal Normal::fromPt(const Pt& a_pt) { return Normal(a_pt); }

inline Normal Normal::fromPtNormalized(const Pt& a_pt) {
  auto normal_to_normalize = Normal(a_pt);
  normal_to_normalize.normalize();
  return normal_to_normalize;
}

inline Pt Normal::toPt(const Normal& a_normal) {
  return Pt(a_normal[0], a_normal[1], a_normal[2]);
}

inline Normal Normal::fromScalarConstant(const double a_constant) {
  return Normal(a_constant);
}

template <class E>
inline Normal::Normal(const Expr<E>& a_expr): normal_m(a_expr) {}

template <class E>
inline Normal::Normal(Expr<E>&& a_expr) : normal_m(std::move(a_expr)) {}

template <class E>
inline Normal& Normal::operator=(const Expr<E>& a_expr) {
  normal_m = a_expr;
  return (*this);
}

template <class E>
inline Normal& Normal::operator=(Expr<E>&& a_expr) {
  normal_m = std::move(a_expr);
  return (*this);
}

inline double& Normal::operator[](const UnsignedIndex_t a_d) {
  assert(a_d < 3);
  return normal_m[a_d];
}

inline const double& Normal::operator[](const UnsignedIndex_t a_d) const {
  assert(a_d < 3);
  return normal_m[a_d];
}

inline constexpr UnsignedIndex_t Normal::size(void) { return 3; }

inline Normal Normal::operator-(void) const {
  return {-normal_m[0], -normal_m[1], -normal_m[2]};
}

inline Normal& Normal::operator+=(const Normal& a_rhs) {
  normal_m[0] += a_rhs[0];
  normal_m[1] += a_rhs[1];
  normal_m[2] += a_rhs[2];
  return (*this);
}

inline Normal& Normal::operator/=(const double a_value) {
  assert(a_value != 0.0);
  for (auto& elem : normal_m) {
    elem /= a_value;
  }
  return (*this);
}

inline Normal& Normal::operator*=(const double a_value) {
  for (auto& elem : normal_m) {
    elem *= a_value;
  }
  return (*this);
}

inline Normal& Normal::operator=(const double a_value) {
  for (auto& elem : normal_m) {
    elem = a_value;
  }
  return (*this);
}

inline bool Normal::operator==(const Normal& a_normal) const {
  return (*this) * a_normal > global_constants::SAME_VEC;
}

inline bool Normal::operator!=(const Normal& a_normal) const {
  return (*this) * a_normal < global_constants::SAME_VEC;
}

inline UnsignedIndex_t Normal::maxMagnitudeComponent(void) const {
  UnsignedIndex_t max = std::fabs((*this)[0]) > std::fabs((*this)[1]) ? 0 : 1;
  max = std::fabs((*this)[max]) > std::fabs((*this)[2]) ? max : 2;
  return max;
}

inline double Normal::calculateMagnitude(void) const {
  return magnitude((*this));
}

inline void Normal::normalize(void) {
  const double inv_magnitude = 1.0 / safelyTiny(this->calculateMagnitude());
  for (auto& elem : normal_m) {
    elem *= inv_magnitude;
  }
}

inline typename Normal::iterator Normal::begin(void) noexcept {
  return normal_m.begin();
}
inline typename Normal::const_iterator Normal::begin(void) const noexcept {
  return this->cbegin();
}
inline typename Normal::const_iterator Normal::cbegin(void) const noexcept {
  return normal_m.cbegin();
}
inline typename Normal::iterator Normal::end(void) noexcept {
  return normal_m.end();
}
inline typename Normal::const_iterator Normal::end(void) const noexcept {
  return this->cend();
}
inline typename Normal::const_iterator Normal::cend(void) const noexcept {
  return normal_m.cend();
}

inline LargeOffsetIndex_t Normal::getSerializedSize(void) const {
  return normal_m.getSerializedSize();
}

inline void Normal::serialize(ByteBuffer* a_buffer) const {
  normal_m.serialize(a_buffer);
}

inline void Normal::unpackSerialized(ByteBuffer* a_buffer) {
  normal_m.unpackSerialized(a_buffer);
}

inline constexpr Normal::Normal(const double* a_normal)
    : normal_m{a_normal[0], a_normal[1], a_normal[2]} {}

inline Normal::Normal(const Pt& a_pt) : normal_m{a_pt[0], a_pt[1], a_pt[2]} {}

inline Normal::Normal(const double a_constant)
    : normal_m{a_constant, a_constant, a_constant} {}

__attribute__((const)) inline double operator*(const Normal& a_normal_0,
                                               const Normal& a_normal_1) {
  return dotProduct(a_normal_0, a_normal_1);
}

__attribute__((const)) inline double operator*(const Pt& a_pt,
                                               const Normal& a_normal) {
  return dotProduct(a_pt, a_normal);
}

__attribute__((const)) inline double operator*(const Normal& a_normal,
                                               const Pt& a_pt) {
  return a_pt * a_normal;
}

__attribute__((const)) inline Normal operator*(const double a_double,
                                               const Normal& a_normal) {
  return Normal(a_normal[0] * a_double, a_normal[1] * a_double,
                a_normal[2] * a_double);
}

__attribute__((const)) inline Normal operator*(const Normal& a_normal,
                                               const double a_double) {
  return a_double * a_normal;
}

__attribute__((const)) inline Normal operator/(const Normal& a_normal,
                                               const double a_double) {
  return Normal(a_normal[0] / a_double, a_normal[1] / a_double,
                a_normal[2] / a_double);
}

inline std::ostream& operator<<(std::ostream& out, const Normal& a_normal) {
  out << std::setprecision(15);
  out << "( " << a_normal[0];
  out << ", " << a_normal[1];
  out << ", " << a_normal[2];
  out << " )";
  return out;
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_GENERAL_NORMAL_TPP_
