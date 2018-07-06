// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_LISTED_VOLUME_MOMENTS_TPP_
#define SRC_MOMENTS_LISTED_VOLUME_MOMENTS_TPP_

namespace IRL {

template <class VolumeMomentsType>
UnsignedIndex_t ListedVolumeMoments<VolumeMomentsType>::size(void) const {
  return list_m.size();
}

template <class VolumeMomentsType>
VolumeMomentsType& ListedVolumeMoments<VolumeMomentsType>::operator[](
    const UnsignedIndex_t a_index) {
  return list_m[a_index];
}

template <class VolumeMomentsType>
const VolumeMomentsType& ListedVolumeMoments<VolumeMomentsType>::operator[](
    const UnsignedIndex_t a_index) const {
  return list_m[a_index];
}

template <class VolumeMomentsType>
ListedVolumeMoments<VolumeMomentsType>& ListedVolumeMoments<VolumeMomentsType>::
operator=(const double a_value) {
  for (auto& moments : (*this)) {
    moments = a_value;
  }
  return (*this);
}

template <class VolumeMomentsType>
void ListedVolumeMoments<VolumeMomentsType>::normalizeByVolume(void) {
  for (auto& volume_moment : (*this)) {
    volume_moment.normalizeByVolume();
  }
}

template <class VolumeMomentsType>
void ListedVolumeMoments<VolumeMomentsType>::multiplyByVolume(void) {
  for (auto& volume_moment : (*this)) {
    volume_moment.multiplyByVolume();
  }
}

template <class VolumeMomentsType>
ListedVolumeMoments<VolumeMomentsType>& ListedVolumeMoments<VolumeMomentsType>::
operator+=(const ListedVolumeMoments& a_rhs) {
  list_m += a_rhs.list_m;
  return (*this);
}

template <class VolumeMomentsType>
ListedVolumeMoments<VolumeMomentsType>& ListedVolumeMoments<VolumeMomentsType>::
operator+=(const VolumeMomentsType& a_rhs) {
  list_m.push_back(a_rhs);
  return (*this);
}

template <class VolumeMomentsType>
void ListedVolumeMoments<VolumeMomentsType>::clear(void) {
  list_m.clear();
}

template <class VolumeMomentsType>
void ListedVolumeMoments<VolumeMomentsType>::erase(
    const UnsignedIndex_t a_index) {
  list_m.erase(a_index);
}

template <class VolumeMomentsType>
typename ListedVolumeMoments<VolumeMomentsType>::iterator
ListedVolumeMoments<VolumeMomentsType>::begin(void) noexcept {
  return list_m.begin();
}
template <class VolumeMomentsType>
typename ListedVolumeMoments<VolumeMomentsType>::const_iterator
ListedVolumeMoments<VolumeMomentsType>::begin(void) const noexcept {
  return this->cbegin();
}
template <class VolumeMomentsType>
typename ListedVolumeMoments<VolumeMomentsType>::const_iterator
ListedVolumeMoments<VolumeMomentsType>::cbegin(void) const noexcept {
  return list_m.cbegin();
}
template <class VolumeMomentsType>
typename ListedVolumeMoments<VolumeMomentsType>::iterator
ListedVolumeMoments<VolumeMomentsType>::end(void) noexcept {
  return list_m.end();
}
template <class VolumeMomentsType>
typename ListedVolumeMoments<VolumeMomentsType>::const_iterator
ListedVolumeMoments<VolumeMomentsType>::end(void) const noexcept {
  return this->cend();
}
template <class VolumeMomentsType>
typename ListedVolumeMoments<VolumeMomentsType>::const_iterator
ListedVolumeMoments<VolumeMomentsType>::cend(void) const noexcept {
  return list_m.cend();
}

template <class VolumeMomentsType>
inline std::ostream& operator<<(
    std::ostream& out,
    const ListedVolumeMoments<VolumeMomentsType>& a_listed_volume_moments) {
  out << std::setprecision(15);
  out << "ListedVolumeMoments contains " << a_listed_volume_moments.size()
      << " members." << std::endl;
  for (const auto& member : a_listed_volume_moments) {
    out << member << std::endl;
  }
  out << std::endl;
  return out;
}
}  // namespace IRL

#endif  // SRC_MOMENTS_LISTED_VOLUME_MOMENTS_TPP_
