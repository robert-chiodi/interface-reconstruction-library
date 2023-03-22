// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2023 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_PT_TPP_
#define IRL_GEOMETRY_GENERAL_PT_TPP_

namespace IRL {

template <class ScalarType>
inline constexpr PtBase<ScalarType>::PtBase(void)
    : loc_m{static_cast<ScalarType>(0), static_cast<ScalarType>(0),
            static_cast<ScalarType>(0)} {}

template <class ScalarType>
inline constexpr PtBase<ScalarType>::PtBase(const UnsignedIndex_t a_x,
                                            const UnsignedIndex_t a_y,
                                            const UnsignedIndex_t a_z)
    : loc_m{static_cast<ScalarType>(a_x), static_cast<ScalarType>(a_y),
            static_cast<ScalarType>(a_z)} {}

template <class ScalarType>
inline constexpr PtBase<ScalarType>::PtBase(const int a_x, const int a_y,
                                            const int a_z)
    : loc_m{static_cast<ScalarType>(a_x), static_cast<ScalarType>(a_y),
            static_cast<ScalarType>(a_z)} {}

template <class ScalarType>
inline constexpr PtBase<ScalarType>::PtBase(const double a_x, const double a_y,
                                            const double a_z)
    : loc_m{a_x, a_y, a_z} {
  static_assert(std::is_same<ScalarType, double>::value,
                "Trying to construct Pt in quad precision with doubles. This "
                "is forbidden to avoid loosing quad precision.");
}

template <class ScalarType>
inline constexpr PtBase<ScalarType>::PtBase(const Quad_t a_x, const Quad_t a_y,
                                            const Quad_t a_z)
    : loc_m{static_cast<ScalarType>(a_x), static_cast<ScalarType>(a_y),
            static_cast<ScalarType>(a_z)} {}

template <class ScalarType>
inline PtBase<ScalarType> PtBase<ScalarType>::fromEdgeIntersection(
    const PtBase<ScalarType>& a_pt_0, const ScalarType a_dist_0,
    const PtBase<ScalarType>& a_pt_1, const ScalarType a_dist_1) {
  const ScalarType mu = a_dist_0 / (a_dist_1 - a_dist_0);
  return a_pt_0 + mu * (a_pt_0 - a_pt_1);
}

template <class ScalarType>
inline PtBase<ScalarType> PtBase<ScalarType>::fromRawDoublePointer(
    const ScalarType* a_loc) {
  return PtBase<ScalarType>(a_loc);
}

template <class ScalarType>
inline PtBase<ScalarType> PtBase<ScalarType>::fromArray(
    const std::array<ScalarType, 3>& a_loc) {
  assert(a_loc.size() == 3);
  return PtBase<ScalarType>(a_loc);
}

template <class ScalarType>
inline constexpr PtBase<ScalarType> PtBase<ScalarType>::fromScalarConstant(
    const ScalarType a_value) {
  return PtBase<ScalarType>(a_value);
}

template <class ScalarType>
inline PtBase<ScalarType> PtBase<ScalarType>::fromVec3(
    const Vec3<ScalarType>& a_vec) {
  return PtBase<ScalarType>(a_vec);
}

template <class ScalarType>
inline PtBase<ScalarType>::operator Vec3<ScalarType>(void) {
  return Vec3<ScalarType>(loc_m[0], loc_m[1], loc_m[2]);
}

template <class ScalarType>
inline PtBase<ScalarType>& PtBase<ScalarType>::operator=(
    const ScalarType a_value) {
  for (auto& loc : loc_m) {
    loc = a_value;
  }
  return (*this);
}

template <class ScalarType>
inline PtBase<ScalarType>& PtBase<ScalarType>::operator*=(
    const ScalarType a_value) {
  for (auto& loc : loc_m) {
    loc *= a_value;
  }
  return (*this);
}

template <class ScalarType>
inline PtBase<ScalarType>& PtBase<ScalarType>::getPt(void) {
  return (*this);
}

template <class ScalarType>
inline const PtBase<ScalarType>& PtBase<ScalarType>::getPt(void) const {
  return (*this);
}

template <class ScalarType>
inline const PtBase<double> PtBase<ScalarType>::toDoublePt(void) const {
  if constexpr (std::is_same<ScalarType, double>::value) {
    return (*this);
  } else {
    return PtBase<double>(static_cast<double>(loc_m[0]),
                          static_cast<double>(loc_m[1]),
                          static_cast<double>(loc_m[2]));
  }
}

template <class ScalarType>
inline const PtBase<Quad_t> PtBase<ScalarType>::toQuadPt(void) const {
  if constexpr (std::is_same<ScalarType, Quad_t>::value) {
    return (*this);
  } else {
    return PtBase<Quad_t>(static_cast<Quad_t>(loc_m[0]),
                          static_cast<Quad_t>(loc_m[1]),
                          static_cast<Quad_t>(loc_m[2]));
  }
}

template <class ScalarType>
inline ScalarType& PtBase<ScalarType>::operator[](const UnsignedIndex_t a_d) {
  assert(a_d < 3);
  return loc_m[a_d];
}

template <class ScalarType>
inline const ScalarType& PtBase<ScalarType>::operator[](
    const UnsignedIndex_t a_d) const {
  assert(a_d < 3);
  return loc_m[a_d];
}

template <class ScalarType>
inline ScalarType& PtBase<ScalarType>::x(void) {
  return (*this)[0];
}

template <class ScalarType>
inline ScalarType& PtBase<ScalarType>::y(void) {
  return (*this)[1];
}

template <class ScalarType>
inline ScalarType& PtBase<ScalarType>::z(void) {
  return (*this)[2];
}

template <class ScalarType>
inline const ScalarType& PtBase<ScalarType>::x(void) const {
  return loc_m[0];
}

template <class ScalarType>
inline const ScalarType& PtBase<ScalarType>::y(void) const {
  return loc_m[1];
}

template <class ScalarType>
inline const ScalarType& PtBase<ScalarType>::z(void) const {
  return loc_m[2];
}

template <>
inline UnsignedIndex_t PtBase<double>::maxComponent(void) const {
  UnsignedIndex_t max = std::fabs((*this)[0]) > std::fabs((*this)[1]) ? 0 : 1;
  max = std::fabs((*this)[max]) > std::fabs((*this)[2]) ? max : 2;
  return max;
}

template <>
inline UnsignedIndex_t PtBase<Quad_t>::maxComponent(void) const {
  UnsignedIndex_t max = fabsq((*this)[0]) > fabsq((*this)[1]) ? 0 : 1;
  max = fabsq((*this)[max]) > fabsq((*this)[2]) ? max : 2;
  return max;
}

template <class ScalarType>
inline PtBase<ScalarType> PtBase<ScalarType>::operator-(void) const {
  return PtBase<ScalarType>(-(*this)[0], -(*this)[1], -(*this)[2]);
}

template <class ScalarType>
inline PtBase<ScalarType>& PtBase<ScalarType>::operator+=(
    const PtBase<ScalarType>& a_rhs) {
  (*this)[0] += a_rhs[0];
  (*this)[1] += a_rhs[1];
  (*this)[2] += a_rhs[2];
  return (*this);
}

template <class ScalarType>
inline PtBase<ScalarType>& PtBase<ScalarType>::operator-=(
    const PtBase<ScalarType>& a_rhs) {
  (*this)[0] -= a_rhs[0];
  (*this)[1] -= a_rhs[1];
  (*this)[2] -= a_rhs[2];
  return (*this);
}

template <class ScalarType>
inline PtBase<ScalarType> PtBase<ScalarType>::operator/(
    const ScalarType a_rhs) {
  auto pt = *this;
  pt /= a_rhs;
  return pt;
}

template <class ScalarType>
inline PtBase<ScalarType>& PtBase<ScalarType>::operator/=(
    const ScalarType a_rhs) {
  (*this)[0] /= a_rhs;
  (*this)[1] /= a_rhs;
  (*this)[2] /= a_rhs;
  return (*this);
}

template <class ScalarType>
inline LargeOffsetIndex_t PtBase<ScalarType>::getSerializedSize(void) const {
  return static_cast<LargeOffsetIndex_t>(loc_m.size() * sizeof(ScalarType));
}

template <class ScalarType>
inline void PtBase<ScalarType>::serialize(ByteBuffer* a_buffer) const {
  a_buffer->pack(loc_m.data(), loc_m.size());
}

template <class ScalarType>
inline void PtBase<ScalarType>::unpackSerialized(ByteBuffer* a_buffer) {
  a_buffer->unpack(loc_m.data(), loc_m.size());
}

template <class ScalarType>
inline constexpr PtBase<ScalarType>::PtBase(const ScalarType* a_loc)
    : loc_m{a_loc[0], a_loc[1], a_loc[2]} {}

template <class ScalarType>
inline constexpr PtBase<ScalarType>::PtBase(
    const std::array<ScalarType, 3>& a_loc)
    : loc_m(a_loc) {}

template <class ScalarType>
inline constexpr PtBase<ScalarType>::PtBase(const ScalarType a_value)
    : loc_m{a_value, a_value, a_value} {}

template <class ScalarType>
inline PtBase<ScalarType>::PtBase(const Vec3<ScalarType>& a_vec)
    : loc_m{a_vec[0], a_vec[1], a_vec[2]} {}

template <class ScalarType>
template <class E>
inline PtBase<ScalarType>::PtBase(const Expr<E>& a_expr) {
  const E& expr(a_expr);
  loc_m[0] = expr[0];
  loc_m[1] = expr[1];
  loc_m[2] = expr[2];
}

template <class ScalarType>
template <class E>
inline PtBase<ScalarType>::PtBase(Expr<E>&& a_expr) {
  const E& expr(a_expr);
  loc_m[0] = expr[0];
  loc_m[1] = expr[1];
  loc_m[2] = expr[2];
}

template <class ScalarType>
template <class E>
inline PtBase<ScalarType>& PtBase<ScalarType>::operator=(
    const Expr<E>& a_expr) {
  const E& expr(a_expr);
  loc_m[0] = expr[0];
  loc_m[1] = expr[1];
  loc_m[2] = expr[2];
  return (*this);
}

template <class ScalarType>
template <class E>
inline PtBase<ScalarType>& PtBase<ScalarType>::operator=(Expr<E>&& a_expr) {
  const E& expr(a_expr);
  loc_m[0] = expr[0];
  loc_m[1] = expr[1];
  loc_m[2] = expr[2];
  return (*this);
}

template <>
inline std::ostream& operator<<(std::ostream& out, const PtBase<double>& a_pt) {
#ifndef NDEBUG
  char* scalar_to_char = new char[30];
  sprintf(scalar_to_char, "%+.20e", a_pt[0]);
  out << "( \033[46m(double) " << scalar_to_char << "\033[0m";
  sprintf(scalar_to_char, "%+.20e", a_pt[1]);
  out << ", \033[46m(double) " << scalar_to_char << "\033[0m";
  sprintf(scalar_to_char, "%+.20e", a_pt[2]);
  out << ", \033[46m(double) " << scalar_to_char << "\033[0m";
  out << " )";
#else
  out << std::setprecision(15);
  out << "( " << a_pt[0];
  out << ", " << a_pt[1];
  out << ", " << a_pt[2];
  out << " )";
#endif
  return out;
}

template <>
inline std::ostream& operator<<(std::ostream& out, const PtBase<Quad_t>& a_pt) {
  out << std::setprecision(15);
  out << "( " << a_pt[0];
  out << ", " << a_pt[1];
  out << ", " << a_pt[2];
  out << " )";
  return out;
}

}  // namespace IRL

#endif  // IRL_GEOMETRY_GENERAL_PT_TPP_
