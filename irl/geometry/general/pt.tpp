// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_PT_TPP_
#define IRL_GEOMETRY_GENERAL_PT_TPP_

namespace IRL {

inline constexpr Pt::Pt(const double a_x, const double a_y, const double a_z)
    : loc_m{a_x, a_y, a_z} {}

inline Pt Pt::fromEdgeIntersection(const Pt& a_pt_0, const double a_dist_0,
                                   const Pt& a_pt_1, const double a_dist_1) {
  const double mu = a_dist_0 / (a_dist_1 - a_dist_0);
  return a_pt_0 + mu * (a_pt_0 - a_pt_1);
}

inline Pt Pt::fromRawDoublePointer(const double* a_loc) { return Pt(a_loc); }

inline Pt Pt::fromArray(const std::array<double, 3>& a_loc) {
  assert(a_loc.size() == 3);
  return Pt(a_loc);
}

inline constexpr Pt Pt::fromScalarConstant(const double a_value) {
  return Pt(a_value);
}

inline Pt Pt::fromVec3(const Vec3& a_vec) { return Pt(a_vec); }

inline Pt::operator Vec3(void) { return Vec3(loc_m[0], loc_m[1], loc_m[2]); }

inline Pt& Pt::operator=(const double a_value) {
  for (auto& loc : loc_m) {
    loc = a_value;
  }
  return (*this);
}

inline Pt& Pt::operator*=(const double a_value) {
  for (auto& loc : loc_m) {
    loc *= a_value;
  }
  return (*this);
}

inline Pt& Pt::getPt(void) { return (*this); }
inline const Pt& Pt::getPt(void) const { return (*this); }

template <class E>
inline Pt::Pt(const Expr<E>& a_expr) {
  const E& expr(a_expr);
  loc_m[0] = expr[0];
  loc_m[1] = expr[1];
  loc_m[2] = expr[2];
}

template <class E>
inline Pt::Pt(Expr<E>&& a_expr) {
  const E& expr(a_expr);
  loc_m[0] = expr[0];
  loc_m[1] = expr[1];
  loc_m[2] = expr[2];
}

template <class E>
inline Pt& Pt::operator=(const Expr<E>& a_expr) {
  const E& expr(a_expr);
  loc_m[0] = expr[0];
  loc_m[1] = expr[1];
  loc_m[2] = expr[2];
  return (*this);
}

template <class E>
inline Pt& Pt::operator=(Expr<E>&& a_expr) {
  const E& expr(a_expr);
  loc_m[0] = expr[0];
  loc_m[1] = expr[1];
  loc_m[2] = expr[2];
  return (*this);
}

inline constexpr UnsignedIndex_t size(const Pt& x) { return 3; }

inline double& Pt::operator[](const UnsignedIndex_t a_d) {
  assert(a_d < 3);
  return loc_m[a_d];
}

inline const double& Pt::operator[](const UnsignedIndex_t a_d) const {
  assert(a_d < 3);
  return loc_m[a_d];
}

inline double& Pt::x(void) { return (*this)[0]; }

inline double& Pt::y(void) { return (*this)[1]; }

inline double& Pt::z(void) { return (*this)[2]; }

inline const double& Pt::x(void) const { return loc_m[0]; }

inline const double& Pt::y(void) const { return loc_m[1]; }

inline const double& Pt::z(void) const { return loc_m[2]; }

inline UnsignedIndex_t Pt::maxComponent(void) const {
  UnsignedIndex_t max = std::fabs((*this)[0]) > std::fabs((*this)[1]) ? 0 : 1;
  max = std::fabs((*this)[max]) > std::fabs((*this)[2]) ? max : 2;
  return max;
}

inline Pt Pt::operator-(void) const {
  return Pt(-(*this)[0], -(*this)[1], -(*this)[2]);
}

inline Pt& Pt::operator+=(const Pt& a_rhs) {
  (*this)[0] += a_rhs[0];
  (*this)[1] += a_rhs[1];
  (*this)[2] += a_rhs[2];
  return (*this);
}

inline Pt& Pt::operator-=(const Pt& a_rhs) {
  (*this)[0] -= a_rhs[0];
  (*this)[1] -= a_rhs[1];
  (*this)[2] -= a_rhs[2];
  return (*this);
}

inline Pt& Pt::operator/(const double a_rhs) {
  (*this)[0] /= a_rhs;
  (*this)[1] /= a_rhs;
  (*this)[2] /= a_rhs;
  return (*this);
}

inline Pt& Pt::operator/=(const double a_rhs) {
  (*this)[0] /= a_rhs;
  (*this)[1] /= a_rhs;
  (*this)[2] /= a_rhs;
  return (*this);
}

inline LargeOffsetIndex_t Pt::getSerializedSize(void) const {
  return static_cast<LargeOffsetIndex_t>(loc_m.size() * sizeof(double));
}

inline void Pt::serialize(ByteBuffer* a_buffer) const {
  a_buffer->pack(loc_m.data(), loc_m.size());
}

inline void Pt::unpackSerialized(ByteBuffer* a_buffer) {
  a_buffer->unpack(loc_m.data(), loc_m.size());
}

inline constexpr Pt::Pt(const double* a_loc)
    : loc_m{a_loc[0], a_loc[1], a_loc[2]} {}

inline constexpr Pt::Pt(const std::array<double, 3>& a_loc) : loc_m(a_loc) {}

inline constexpr Pt::Pt(const double a_value)
    : loc_m{a_value, a_value, a_value} {}

inline Pt::Pt(const Vec3& a_vec) : loc_m{a_vec[0], a_vec[1], a_vec[2]} {}

inline std::ostream& operator<<(std::ostream& out, const Pt& a_pt) {
  out << std::setprecision(15);
  out << "( " << a_pt[0];
  out << ", " << a_pt[1];
  out << ", " << a_pt[2];
  out << " )";
  return out;
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_GENERAL_PT_TPP_
