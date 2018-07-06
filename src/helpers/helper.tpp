// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_HELPERS_HELPER_TPP_
#define SRC_HELPERS_HELPER_TPP_

namespace IRL {

inline double safelyTiny(const double a_value) {
  return std::copysign(std::max(std::fabs(a_value), DBL_MIN), a_value);
}

inline double safelyEpsilon(const double a_value) {
  return std::copysign(std::max(std::fabs(a_value), DBL_EPSILON), a_value);
}

inline double safelySmall(const double a_value, const double a_small_value) {
  assert(a_small_value > 0.0);
  return std::copysign(std::max(std::fabs(a_value), a_small_value), a_value);
}

inline double clipBetween(const double a_smallest_value, const double a_value,
                          const double a_largest_value) {
  return std::max(a_smallest_value, std::min(a_value, a_largest_value));
}

inline bool wantPurelyInternal(const double a_internal_fraction) {
  return a_internal_fraction > global_constants::VF_HIGH;
}

inline bool wantPurelyExternal(const double a_internal_fraction) {
  return a_internal_fraction < global_constants::VF_LOW;
}

inline void sort3Ascending(double* a_items) {
  double tmp;
  if (a_items[1] < a_items[0]) {
    tmp = a_items[1];
    a_items[1] = a_items[0];
    a_items[0] = tmp;
  }
  if (a_items[2] < a_items[1]) {
    tmp = a_items[2];
    a_items[2] = a_items[1];
    a_items[1] = tmp;
  }
  if (a_items[1] < a_items[0]) {
    tmp = a_items[1];
    a_items[1] = a_items[0];
    a_items[0] = tmp;
  }
}

inline void sort3Descending(double* a_items) {
  double tmp;
  if (a_items[1] > a_items[0]) {
    tmp = a_items[1];
    a_items[1] = a_items[0];
    a_items[0] = tmp;
  }
  if (a_items[2] > a_items[1]) {
    tmp = a_items[2];
    a_items[2] = a_items[1];
    a_items[1] = tmp;
  }
  if (a_items[1] > a_items[0]) {
    tmp = a_items[1];
    a_items[1] = a_items[0];
    a_items[0] = tmp;
  }
}

template <class IndexType>
inline void sort3AscendingTracked(double* a_items,
                                  IndexType* a_original_index) {
  double tmp;
  IndexType tmpi;
  if (a_items[1] < a_items[0]) {
    tmp = a_items[1];
    a_items[1] = a_items[0];
    a_items[0] = tmp;
    tmpi = a_original_index[1];
    a_original_index[1] = a_original_index[0];
    a_original_index[0] = tmpi;
  }
  if (a_items[2] < a_items[1]) {
    tmp = a_items[2];
    a_items[2] = a_items[1];
    a_items[1] = tmp;
    tmpi = a_original_index[2];
    a_original_index[2] = a_original_index[1];
    a_original_index[1] = tmpi;
  }
  if (a_items[1] < a_items[0]) {
    tmp = a_items[1];
    a_items[1] = a_items[0];
    a_items[0] = tmp;
    tmpi = a_original_index[1];
    a_original_index[1] = a_original_index[0];
    a_original_index[0] = tmpi;
  }
}

template <class IndexType>
inline void sort3DescendingTracked(double* a_items,
                                   IndexType* a_original_index) {
  double tmp;
  IndexType tmpi;
  if (a_items[1] > a_items[0]) {
    tmp = a_items[1];
    a_items[1] = a_items[0];
    a_items[0] = tmp;
    tmpi = a_original_index[1];
    a_original_index[1] = a_original_index[0];
    a_original_index[0] = tmpi;
  }
  if (a_items[2] > a_items[1]) {
    tmp = a_items[2];
    a_items[2] = a_items[1];
    a_items[1] = tmp;
    tmpi = a_original_index[2];
    a_original_index[2] = a_original_index[1];
    a_original_index[1] = tmpi;
  }
  if (a_items[1] > a_items[0]) {
    tmp = a_items[1];
    a_items[1] = a_items[0];
    a_items[0] = tmp;
    tmpi = a_original_index[1];
    a_original_index[1] = a_original_index[0];
    a_original_index[0] = tmpi;
  }
}

template <class CarriedType, class DictatingType>
inline void sortAscendingBasedOnOtherArray(CarriedType* a_carried_array,
                                           DictatingType* a_dictating_array) {
  assert(a_carried_array->size() == a_dictating_array->size());
  for (UnsignedIndex_t i = 1; i < a_dictating_array->size(); ++i) {
    for (UnsignedIndex_t j = i;
         j != 0 && (*a_dictating_array)[j - 1] > (*a_dictating_array)[j]; --j) {
      std::swap((*a_carried_array)[j], (*a_carried_array)[j - 1]);
      std::swap((*a_dictating_array)[j], (*a_dictating_array)[j - 1]);
    }
  }
}

template <class CarriedType, class DictatingType>
inline void sortDescendingBasedOnOtherArray(CarriedType* a_carried_array,
                                            DictatingType* a_dictating_array) {
  assert(a_carried_array->size() == a_dictating_array->size());
  for (UnsignedIndex_t i = 1; i < a_dictating_array->size(); ++i) {
    for (UnsignedIndex_t j = i;
         j != 0 && (*a_dictating_array)[j - 1] < (*a_dictating_array)[j]; --j) {
      std::swap((*a_carried_array)[j], (*a_carried_array)[j - 1]);
      std::swap((*a_dictating_array)[j], (*a_dictating_array)[j - 1]);
    }
  }
}

}  // namespace IRL

#endif  // SRC_HELPERS_HELPER_TPP_
