// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_HELPERS_MYMATH_H_
#define IRL_HELPERS_MYMATH_H_

#include <algorithm>
#include <cmath>
#include <ostream>
#include <string>

#include "quadmath.h"

#include "irl/helpers/helper.h"
#include "irl/parameters/defined_types.h"

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
inline typename DataType::value_type magnitude(const DataType& a_vector);

/// \brief Calculate the squared magnitude between for a 3 element vector.
///
/// Template requirements for DataType:
/// - Overload `operator[]` : Method that supplies
/// access to the underlying variable that cross products
/// will be taken against.
template <class DataType>
inline typename DataType::value_type squaredMagnitude(const DataType& a_vector);

/// \brief Take cross product of two 3-element vectors
///
/// Template requirements for DataType:
/// - Overload `operator[]` : Method that supplies
/// access to the underlying variable to use for calculating magnitude.
// template <class DataType>
// inline vector_cross_product crossProduct(const DataType& a_vector_0,
//                              const DataType& a_vector_1);

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
inline typename T1::value_type dotProduct(const T1& a_vector_0,
                                          const T2& a_vector_1);

/// \brief Scalar triple product of 3, 3 element vectors
///
/// Template requirements for DataType:
/// - Overload `operator[]` : Method that supplies
/// access to the underlying variable that cross products
/// will be taken against.
template <class DataType>
inline typename DataType::value_type scalarTripleProduct(
    const DataType& a_vector_0, const DataType& a_vector_1,
    const DataType& a_vector_2);

////////// Overloading math functions to allow for quad precision
inline std::ostream& operator<<(std::ostream& out, const Quad_t a_scalar);

template <class ScalarType>
inline ScalarType machine_epsilon(void);

template <class ScalarType>
inline ScalarType machine_pi(void);

template <class ScalarType>
inline ScalarType abs(const ScalarType a_scalar);

template <class ScalarType>
inline ScalarType fabs(const ScalarType a_scalar);

template <class ScalarType>
inline ScalarType sqrt(const ScalarType a_scalar);

template <class ScalarType, class PowerType>
inline ScalarType pow(const ScalarType a_scalar, const PowerType a_power);

template <class ScalarType>
inline ScalarType log(const ScalarType a_scalar);

template <class ScalarType>
inline ScalarType atan(const ScalarType a_scalar);

template <class ScalarType>
inline ScalarType atan2(const ScalarType a_scalar_y,
                        const ScalarType a_scalar_x);

template <class ScalarType>
inline ScalarType atanh(const ScalarType a_scalar);

template <class ScalarType>
inline ScalarType cos(const ScalarType a_scalar);

template <class ScalarType>
inline ScalarType sin(const ScalarType a_scalar);

template <class ScalarType>
inline ScalarType copysign(const ScalarType a_scalar, const ScalarType a_sign);

template <class ScalarType>
inline ScalarType minimum(const ScalarType a_scalar1,
                          const ScalarType a_scalar2);

template <class ScalarType>
inline ScalarType maximum(const ScalarType a_scalar1,
                          const ScalarType a_scalar2);

}  // namespace IRL

#include "irl/helpers/mymath.tpp"

#endif  // IRL_HELPERS_MYMATH_H_
