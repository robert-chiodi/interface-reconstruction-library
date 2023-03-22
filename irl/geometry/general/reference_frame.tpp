// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_REFERENCE_FRAME_TPP_
#define IRL_GEOMETRY_GENERAL_REFERENCE_FRAME_TPP_

namespace IRL {

template <class ScalarType>
inline ReferenceFrameBase<ScalarType>::ReferenceFrameBase(
    const NormalBase<ScalarType>& a_axis_0,
    const NormalBase<ScalarType>& a_axis_1,
    const NormalBase<ScalarType>& a_axis_2)
    : axis_m{a_axis_0, a_axis_1, a_axis_2} {}

template <class ScalarType>
inline NormalBase<ScalarType>& ReferenceFrameBase<ScalarType>::operator[](
    const UnsignedIndex_t a_axis) {
  assert(a_axis < 3);
  return axis_m[a_axis];
}

template <class ScalarType>
inline const NormalBase<ScalarType>& ReferenceFrameBase<ScalarType>::operator[](
    const UnsignedIndex_t a_axis) const {
  assert(a_axis < 3);
  return axis_m[a_axis];
}

template <class ScalarType>
inline bool ReferenceFrameBase<ScalarType>::isOrthonormalBasis(void) const {
  static constexpr ScalarType TOLERANCE =
      static_cast<ScalarType>(10.0 * DBL_EPSILON);
  if (fabs(axis_m[0] * axis_m[1]) > TOLERANCE) {
    return false;
  }
  if (fabs(axis_m[0] * axis_m[2]) > TOLERANCE) {
    return false;
  }
  if (fabs(axis_m[1] * axis_m[2]) > TOLERANCE) {
    return false;
  }
  return true;
}

template <class ScalarType>
inline typename ReferenceFrameBase<ScalarType>::iterator
ReferenceFrameBase<ScalarType>::begin(void) noexcept {
  return axis_m.begin();
}

template <class ScalarType>
inline typename ReferenceFrameBase<ScalarType>::const_iterator
ReferenceFrameBase<ScalarType>::begin(void) const noexcept {
  return this->cbegin();
}

template <class ScalarType>
inline typename ReferenceFrameBase<ScalarType>::const_iterator
ReferenceFrameBase<ScalarType>::cbegin(void) const noexcept {
  return axis_m.cbegin();
}

template <class ScalarType>
inline typename ReferenceFrameBase<ScalarType>::iterator
ReferenceFrameBase<ScalarType>::end(void) noexcept {
  return axis_m.end();
}

template <class ScalarType>
inline typename ReferenceFrameBase<ScalarType>::const_iterator
ReferenceFrameBase<ScalarType>::end(void) const noexcept {
  return this->cend();
}

template <class ScalarType>
inline typename ReferenceFrameBase<ScalarType>::const_iterator
ReferenceFrameBase<ScalarType>::cend(void) const noexcept {
  return axis_m.cend();
}

}  // namespace IRL

#endif  // IRL_GEOMETRY_GENERAL_REFERENCE_FRAME_TPP_
