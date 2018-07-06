// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_ACCUMULATED_VOLUME_MOMENTS_TPP_
#define SRC_MOMENTS_ACCUMULATED_VOLUME_MOMENTS_TPP_

namespace IRL {

template <class VolumeMomentsType>
UnsignedIndex_t AccumulatedVolumeMoments<VolumeMomentsType>::size(void) const {
  return accumulated_moments_m.size();
}

template <class VolumeMomentsType>
VolumeMomentsType& AccumulatedVolumeMoments<VolumeMomentsType>::operator[](
    const UnsignedIndex_t a_index) {
  return accumulated_moments_m[a_index];
}

template <class VolumeMomentsType>
const VolumeMomentsType& AccumulatedVolumeMoments<VolumeMomentsType>::
operator[](const UnsignedIndex_t a_index) const {
  return accumulated_moments_m[a_index];
}

template <class VolumeMomentsType>
AccumulatedVolumeMoments<VolumeMomentsType>&
AccumulatedVolumeMoments<VolumeMomentsType>::operator=(const double a_value) {
  for (auto& moments : (*this)) {
    moments = a_value;
  }
  return (*this);
}

template <class VolumeMomentsType>
void AccumulatedVolumeMoments<VolumeMomentsType>::normalizeByVolume(void) {
  for (auto& volume_moment : (*this)) {
    volume_moment.normalizeByVolume();
  }
}

template <class VolumeMomentsType>
void AccumulatedVolumeMoments<VolumeMomentsType>::multiplyByVolume(void) {
  for (auto& volume_moment : (*this)) {
    volume_moment.multiplyByVolume();
  }
}

template <class VolumeMomentsType>
AccumulatedVolumeMoments<VolumeMomentsType>&
AccumulatedVolumeMoments<VolumeMomentsType>::operator+=(
    const AccumulatedVolumeMoments& a_rhs) {
  accumulated_moments_m += a_rhs.accumulated_moments_m;
  return (*this);
}

template <class VolumeMomentsType>
void AccumulatedVolumeMoments<VolumeMomentsType>::clear(void) {
  accumulated_moments_m.clear();
}

template <class VolumeMomentsType>
typename AccumulatedVolumeMoments<VolumeMomentsType>::iterator
AccumulatedVolumeMoments<VolumeMomentsType>::begin(void) noexcept {
  return accumulated_moments_m.begin();
}
template <class VolumeMomentsType>
typename AccumulatedVolumeMoments<VolumeMomentsType>::const_iterator
AccumulatedVolumeMoments<VolumeMomentsType>::begin(void) const noexcept {
  return this->cbegin();
}
template <class VolumeMomentsType>
typename AccumulatedVolumeMoments<VolumeMomentsType>::const_iterator
AccumulatedVolumeMoments<VolumeMomentsType>::cbegin(void) const noexcept {
  return accumulated_moments_m.cbegin();
}
template <class VolumeMomentsType>
typename AccumulatedVolumeMoments<VolumeMomentsType>::iterator
AccumulatedVolumeMoments<VolumeMomentsType>::end(void) noexcept {
  return accumulated_moments_m.end();
}
template <class VolumeMomentsType>
typename AccumulatedVolumeMoments<VolumeMomentsType>::const_iterator
AccumulatedVolumeMoments<VolumeMomentsType>::end(void) const noexcept {
  return this->cend();
}
template <class VolumeMomentsType>
typename AccumulatedVolumeMoments<VolumeMomentsType>::const_iterator
AccumulatedVolumeMoments<VolumeMomentsType>::cend(void) const noexcept {
  return accumulated_moments_m.cend();
}

template <class VolumeMomentsType>
inline AccumulatedVolumeMoments<VolumeMomentsType> operator*(
    const AccumulatedVolumeMoments<VolumeMomentsType>& a_list,
    const double a_multiplier) {
  auto list_to_return = a_list;
  for (auto& member : list_to_return) {
    member *= a_multiplier;
  }
  return list_to_return;
}
template <class VolumeMomentsType>
inline AccumulatedVolumeMoments<VolumeMomentsType> operator*(
    const double a_multiplier,
    const AccumulatedVolumeMoments<VolumeMomentsType>& a_list) {
  return a_list * a_multiplier;
}

}  // namespace IRL

#endif  // SRC_MOMENTS_ACCUMULATED_VOLUME_MOMENTS_TPP_
