// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_HELPERS_MYMATH_H_
#define SRC_HELPERS_MYMATH_H_

#include <algorithm>
#include <cmath>

#include "src/helpers/helper.h"

namespace IRL {

/// \file mymath.h
///
/// This file contains functions that deal with
/// math like operations like conversions of
/// angle types, dot products, cross products,
/// etc.
///
/// First, the function declarations are given.
/// Afterwards, the inlined function definitions are given.
/// Then, templated functions are given.

/// \brief Convert degrees to radians.
inline constexpr double deg2Rad(const double a_degree);

/// \brief Convert radians to degrees.
inline constexpr double rad2Deg(const double a_radian);

/// \brief Normalize a given angle (in radians) to be between 0 and 2\f$pi\f$.
inline double angleNormalize(const double a_radian);

/// \brief Normalize a given angle (in radians) to be between -2\f$pi\f$ and
/// 2\f$pi\f$.
inline double signedAngleNormalize(const double a_radian);

/// \brief Calculate the magnitude between two 3 element vectors
///
/// Template requirements for DataType:
/// - Overload `operator[]` : Method that supplies
/// access to the underlying variable to use for calculating magnitude.
template <class DataType>
inline double magnitude(const DataType& a_vector_0, const DataType& a_vector_1);

/// \brief Calculate the squared magnitude between for a 3 element vector.
///
/// Template requirements for DataType:
/// - Overload `operator[]` : Method that supplies
/// access to the underlying variable that cross products
/// will be taken against.
template <class DataType>
inline double squaredMagnitude(const DataType& a_vector);

/// \brief Take cross product of two 3-element vectors
///
/// Template requirements for DataType:
/// - Overload `operator[]` : Method that supplies
/// access to the underlying variable to use for calculating magnitude.
//template <class DataType>
//inline vector_cross_product crossProduct(const DataType& a_vector_0,
//                             const DataType& a_vector_1);

/// \brief Cross product of two 3-element vectors then normalized
///
/// Template requirements for DataType:
/// - Overload `operator[]` : Method that supplies
/// access to the underlying variable that cross products
/// will be taken against.
/// - Overload `operator/=` : Method that divides
/// each of its three components by a double to
/// normalize the vector.
template <class DataType>
inline DataType crossProductNormalized(const DataType& a_vector_0,
                                       const DataType& a_vector_1);

/// \brief Dot product between two 3 element vectors.
///
/// Template requirements for DataType:
/// - Overload `operator[]` : Method that supplies
/// access to the underlying variable that cross products
/// will be taken against.
template <class T1, class T2>
inline double dotProduct(const T1& a_vector_0,
                         const T2& a_vector_1);

/// \brief Scalar triple product of 3, 3 element vectors
///
/// Template requirements for DataType:
/// - Overload `operator[]` : Method that supplies
/// access to the underlying variable that cross products
/// will be taken against.
template <class DataType>
inline double scalarTripleProduct(const DataType& a_vector_0,
                                  const DataType& a_vector_1,
                                  const DataType& a_vector_2);

}  // namespace IRL

#include "src/helpers/mymath.tpp"

#endif  // SRC_HELPERS_MYMATH_H_
