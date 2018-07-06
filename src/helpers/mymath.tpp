// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_HELPERS_MYMATH_TPP_
#define SRC_HELPERS_MYMATH_TPP_

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
inline double magnitude(const DataType& a_vector) {
  return std::sqrt(a_vector[0] * a_vector[0] + a_vector[1] * a_vector[1] +
                   a_vector[2] * a_vector[2]);
}

template <class DataType>
inline double squaredMagnitude(const DataType& a_vector) {
  return a_vector[0] * a_vector[0] + a_vector[1] * a_vector[1] +
         a_vector[2] * a_vector[2];
}

//template <class DataType>
//inline DataType crossProduct(const DataType& a_vector_0,
//                             const DataType& a_vector_1) {
//  return DataType(
//      a_vector_0[1] * a_vector_1[2] - a_vector_0[2] * a_vector_1[1],
//      a_vector_0[2] * a_vector_1[0] - a_vector_0[0] * a_vector_1[2],
//      a_vector_0[0] * a_vector_1[1] - a_vector_0[1] * a_vector_1[0]);
//}

template <class DataType>
inline DataType crossProductNormalized(const DataType& a_vector_0,
                                       const DataType& a_vector_1) {
  DataType object_to_return = crossProduct(a_vector_0, a_vector_1);
  object_to_return /= safelyTiny(magnitude(object_to_return));
  return object_to_return;
}

template <class T1, class T2>
inline double dotProduct(const T1& a_vector_0,
                         const T2& a_vector_1) {
  return a_vector_0[0] * a_vector_1[0] + a_vector_0[1] * a_vector_1[1] +
         a_vector_0[2] * a_vector_1[2];
}

template <class DataType>
inline double scalarTripleProduct(const DataType& a_vector_0,
                                  const DataType& a_vector_1,
                                  const DataType& a_vector_2) {
  return dotProduct(a_vector_0, crossProduct(a_vector_1, a_vector_2));
}

}  // namespace IRL

#endif  // SRC_HELPERS_MYMATH_TPP_
