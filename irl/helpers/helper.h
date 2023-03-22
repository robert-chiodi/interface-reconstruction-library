// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_HELPERS_HELPER_H_
#define IRL_HELPERS_HELPER_H_

#include <float.h>
#include <quadmath.h>
#include <type_traits>
#include <utility>

#include <algorithm>
#include <cassert>
#include <cmath>

#include "irl/helpers/SFINAE_boiler_plate.h"
#include "irl/helpers/expression_templates.h"
#include "irl/parameters/constants.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

/// \file helper.h
///
/// Helper functions that perform small
/// but necessary functions that make coding
/// simpler and prevents code bloat.
///
/// First, the function declarations are given.
/// Afterwards, the inlined function definitions are given.

/// \brief Takes max between abs(`a_value`) and DBL_MIN while preserving sign.
inline double safelyTiny(const double a_value);
inline Quad_t safelyTiny(const Quad_t a_value);

/// \brief Takes max between abs(`a_value`) and DBL_EPSILON while preserving
/// sign.
inline double safelyEpsilon(const double a_value);
inline Quad_t safelyEpsilon(const Quad_t a_value);

/// \brief Takes max between abs(`a_value`) and `the_smallest_value` while
/// preserving sign.
inline double safelySmall(const double a_value, const double a_small_value);

/// \brief Clip value to lay between the min and max values
inline double clipBetween(const double a_smallest_value, const double a_value,
                          const double a_largest_value);

/// \brief Returns whether the liquid volume fraction indicates full liquid
inline bool wantPurelyInternal(const double a_internal_fraction);

/// \brief Returns whether the liquid volume fraction indicates full gas
inline bool wantPurelyExternal(const double a_internal_fraction);

/// \brief Sort 3 doubles in an array into ascending order via insertion.
inline void sort3Ascending(double* a_items);

/// \brief Sort 3 doubles in an array into descending order via insertion.
inline void sort3Descending(double* a_items);

/// \brief Sort 3 doubles in an array into ascending order via insertion.
/// The original location of the item is also tracked to allow reversing.
template <class IndexType>
inline void sort3AscendingTracked(double* a_items, IndexType* a_original_index);

/// \brief Sort 3 doubles in an array into descending order via insertion.
/// The original location of the item is also tracked to allow reversing.
template <class IndexType>
inline void sort3DescendingTracked(double* a_items,
                                   IndexType* a_original_index);

/// \brief This function sorts `a_carried_array` based on the
/// elements in `a_dictating_array` using insertion sort.
///
/// Template requirements for CarriedType:
/// - Access via overloaded operator[].
/// - Ability to be used with std::swap().
///
/// Template requirements for DictatingType:
/// - Access via overloaded operator[].
/// - Ability to be used with std::swap().
/// - Comparison operator >(int i1, int i2) that returns whether one
/// the element at i1 is greater than the element at i2.
///
/// @param[inout] a_carried_array Array that will be sorted into ascending
/// order according to the elements in `a_dictating_array`
/// @param[inout] a_dictating_array Array that will be sorted and
/// drive the sorting in `a_carried_array`
/// @param[in] a_length_of_array Length of both arrays being sorted.
template <class CarriedType, class DictatingType>
inline void sortAscendingBasedOnOtherArray(CarriedType* a_carried_array,
                                           DictatingType* a_dictating_array);

/// \brief This function sorts `a_carried_array` based on the
/// elements in `a_dictating_array` using insertion sort.
///
/// Template requirements for CarriedType:
/// - Access via overloaded operator[].
/// - Ability to be used with std::swap().
///
/// Template requirements for DictatingType:
/// - Access via overloaded operator[].
/// - Ability to be used with std::swap().
/// - Comparison operator <(int i1, int i2) that returns whether one
/// the element at i1 is less than the element at i2.
///
/// @param[inout] a_carried_array Array that will be sorted into descending
/// order according to the elements in `a_dictating_array`
/// @param[inout] a_dictating_array Array that will be sorted and
/// drive the sorting in `a_carried_array`
/// @param[in] a_length_of_array Length of both arrays being sorted.
template <class CarriedType, class DictatingType>
inline void sortDescendingBasedOnOtherArray(CarriedType* a_carried_array,
                                            DictatingType* a_dictating_array);

}  // namespace IRL

#include "irl/helpers/helper.tpp"

#endif  // IRL_HELPERS_HELPER_H_
