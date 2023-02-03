// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_HELPERS_MYMATH_TPP_
#define IRL_HELPERS_MYMATH_TPP_

#ifdef __x86_64__
#include <immintrin.h>
#else
#include "sse2neon.h"
#endif

namespace IRL {

inline constexpr double deg2Rad(const double a_degree) {
  return a_degree * M_PI / 180.0;
}

inline constexpr double rad2Deg(const double a_radian) {
  return a_radian * 180.0 / M_PI;
}

inline double angleNormalize(const double a_radian) {
  return std::max(0.0,
                  a_radian - (std::floor(0.5 * a_radian / M_PI)) * 2.0 * M_PI);
}
inline double signedAngleNormalize(const double a_radian) {
  return a_radian < 0.0 ? -(2.0 * M_PI - angleNormalize(a_radian))
                        : angleNormalize(a_radian);
}

template <class DataType>
inline typename DataType::value_type magnitude(const DataType& a_vector) {
  return sqrt(a_vector[0] * a_vector[0] + a_vector[1] * a_vector[1] +
              a_vector[2] * a_vector[2]);
}

template <class DataType>
inline typename DataType::value_type squaredMagnitude(
    const DataType& a_vector) {
  return a_vector[0] * a_vector[0] + a_vector[1] * a_vector[1] +
         a_vector[2] * a_vector[2];
}

// template <class DataType>
// inline DataType crossProduct(const DataType& a_vector_0,
//                              const DataType& a_vector_1) {
//   return DataType(
//       a_vector_0[1] * a_vector_1[2] - a_vector_0[2] * a_vector_1[1],
//       a_vector_0[2] * a_vector_1[0] - a_vector_0[0] * a_vector_1[2],
//       a_vector_0[0] * a_vector_1[1] - a_vector_0[1] * a_vector_1[0]);
// }

template <class DataType>
inline DataType crossProductNormalized(const DataType& a_vector_0,
                                       const DataType& a_vector_1) {
  DataType object_to_return = crossProduct(a_vector_0, a_vector_1);
  object_to_return /= safelyTiny(magnitude(object_to_return));
  return object_to_return;
}

template <class T1, class T2>
inline typename T1::value_type dotProduct(const T1& a_vector_0,
                                          const T2& a_vector_1) {
  static_assert(
      std::is_same<typename T1::value_type, typename T2::value_type>::value,
      "Trying to dot two vectors of incompatible types (double and quad).");
  return a_vector_0[0] * a_vector_1[0] + a_vector_0[1] * a_vector_1[1] +
         a_vector_0[2] * a_vector_1[2];
}

template <class DataType>
inline typename DataType::value_type scalarTripleProduct(
    const DataType& a_vector_0, const DataType& a_vector_1,
    const DataType& a_vector_2) {
  return dotProduct(a_vector_0, crossProduct(a_vector_1, a_vector_2));
}

// Overloading math functions to allow for quad precision
inline std::ostream& operator<<(std::ostream& out, const Quad_t a_scalar) {
#ifndef NDEBUG
  char* scalar_to_char = new char[30];
  quadmath_snprintf(scalar_to_char, 30, "%+.20Qe", a_scalar);
  out << "\033[44m(__float128) " << scalar_to_char << "\033[0m";
#else
  out << static_cast<double>(a_scalar);
#endif
  return out;
}

template <>
inline bool isnan(double a_scalar) {
  return std::isnan(a_scalar);
}

template <>
inline bool isnan(Quad_t a_scalar) {
  return isnanq(a_scalar) == 1;
}

template <>
inline double machine_epsilon(void) {
  return DBL_EPSILON;
}

template <>
inline Quad_t machine_epsilon(void) {
  return 1.0e5q * FLT128_EPSILON;
}

template <>
inline double machine_pi(void) {
  return M_PI;
}

template <>
inline Quad_t machine_pi(void) {
  return M_PIq;
}

template <>
inline double abs(const double a_scalar) {
  return std::fabs(a_scalar);
}

template <>
inline Quad_t abs(const Quad_t a_scalar) {
  return fabsq(a_scalar);
}

template <>
inline double fabs(const double a_scalar) {
  return std::fabs(a_scalar);
}

template <>
inline Quad_t fabs(const Quad_t a_scalar) {
  return fabsq(a_scalar);
}

template <>
inline double sqrt(const double a_scalar) {
  return std::sqrt(a_scalar);
}

template <>
inline Quad_t sqrt(const Quad_t a_scalar) {
  return sqrtq(a_scalar);
}

template <>
inline double approxinvsqrt(const double a_scalar) {
  __m128 temp = _mm_set_ss(static_cast<float>(a_scalar));
  temp = _mm_rsqrt_ss(temp);
  return static_cast<double>(_mm_cvtss_f32(temp));
}

template <>
inline Quad_t approxinvsqrt(const Quad_t a_scalar) {
  return 1.0q / sqrtq(a_scalar);
}

template <>
inline double invsqrt(const double a_scalar) {
  return 1.0 / std::sqrt(a_scalar);
}

template <>
inline Quad_t invsqrt(const Quad_t a_scalar) {
  return 1.0q / sqrtq(a_scalar);
}

template <>
inline double pow(const double a_scalar, const double a_power) {
  return std::pow(a_scalar, a_power);
}

template <>
inline Quad_t pow(const Quad_t a_scalar, const int a_power) {
  return powq(a_scalar, static_cast<Quad_t>(a_power));
}

template <>
inline Quad_t pow(const Quad_t a_scalar, const double a_power) {
  return powq(a_scalar, static_cast<Quad_t>(a_power));
}

template <>
inline Quad_t pow(const Quad_t a_scalar, const Quad_t a_power) {
  return powq(a_scalar, a_power);
}

template <>
inline double log(const double a_scalar) {
  return std::log(a_scalar);
}

template <>
inline Quad_t log(const Quad_t a_scalar) {
  return logq(a_scalar);
}

template <>
inline double atan2(const double a_scalar_y, const double a_scalar_x) {
  return std::atan2(a_scalar_y, a_scalar_x);
}

template <>
inline Quad_t atan2(const Quad_t a_scalar_y, const Quad_t a_scalar_x) {
  return atan2q(a_scalar_y, a_scalar_x);
}

template <>
inline double atan(const double a_scalar) {
  return std::atan(a_scalar);
}

template <>
inline Quad_t atan(const Quad_t a_scalar) {
  return atanq(a_scalar);
}

template <>
inline double atanh(const double a_scalar) {
  return std::atanh(a_scalar);
}

template <>
inline Quad_t atanh(const Quad_t a_scalar) {
  return atanhq(a_scalar);
}

template <>
inline double cos(const double a_scalar) {
  return std::cos(a_scalar);
}

template <>
inline Quad_t cos(const Quad_t a_scalar) {
  return cosq(a_scalar);
}

template <>
inline double sin(const double a_scalar) {
  return std::sin(a_scalar);
}

template <>
inline Quad_t sin(const Quad_t a_scalar) {
  return sinq(a_scalar);
}

template <>
inline double copysign(const double a_scalar, const double a_sign) {
  return std::copysign(a_scalar, a_sign);
}

template <>
inline Quad_t copysign(const Quad_t a_scalar, const Quad_t a_sign) {
  if (a_sign >= 0.0q) {
    return fabs(a_scalar);
  } else {
    return -fabs(a_scalar);
  }
}

template <>
inline double minimum(const double a_scalar1, const double a_scalar2) {
  return std::min(a_scalar1, a_scalar2);
}

template <>
inline Quad_t minimum(const Quad_t a_scalar1, const Quad_t a_scalar2) {
  if (a_scalar1 < a_scalar2) {
    return a_scalar1;
  } else {
    return a_scalar2;
  }
}

template <>
inline double maximum(const double a_scalar1, const double a_scalar2) {
  return std::max(a_scalar1, a_scalar2);
}

template <>
inline Quad_t maximum(const Quad_t a_scalar1, const Quad_t a_scalar2) {
  if (a_scalar1 >= a_scalar2) {
    return a_scalar1;
  } else {
    return a_scalar2;
  }
}

}  // namespace IRL

#endif  // IRL_HELPERS_MYMATH_TPP_
