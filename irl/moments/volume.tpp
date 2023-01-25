// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_VOLUME_TPP_
#define IRL_MOMENTS_VOLUME_TPP_

namespace IRL {

template <class ScalarType>
inline VolumeBase<ScalarType>::VolumeBase(void)
    : volume_m{static_cast<ScalarType>(0)} {}

template <class ScalarType>
inline constexpr VolumeBase<ScalarType>::VolumeBase(const ScalarType a_value)
    : volume_m{a_value} {}

template <class ScalarType>
inline constexpr VolumeBase<ScalarType>::VolumeBase(
    const VolumeBase<double>& a_value)
    : volume_m{static_cast<ScalarType>(a_value.volume())} {}

template <class ScalarType>
inline constexpr VolumeBase<ScalarType>::VolumeBase(
    const VolumeBase<Quad_t>& a_value)
    : volume_m{static_cast<ScalarType>(a_value.volume())} {}

template <class ScalarType>
inline VolumeBase<ScalarType> VolumeBase<ScalarType>::fromScalarConstant(
    const ScalarType a_value) {
  return VolumeBase<ScalarType>(a_value);
}

template <class ScalarType>
template <class GeometryType>
inline VolumeBase<ScalarType> VolumeBase<ScalarType>::calculateMoments(
    GeometryType* a_geometry) {
  return a_geometry->calculateVolume();
}

template <class ScalarType>
inline void VolumeBase<ScalarType>::multiplyByVolume(void) {}

template <class ScalarType>
inline void VolumeBase<ScalarType>::normalizeByVolume(void) {}

template <class ScalarType>
inline ScalarType& VolumeBase<ScalarType>::volume(void) {
  return volume_m;
}

template <class ScalarType>
inline const ScalarType& VolumeBase<ScalarType>::volume(void) const {
  return volume_m;
}

template <class ScalarType>
inline VolumeBase<ScalarType>& VolumeBase<ScalarType>::operator+=(
    const VolumeBase<ScalarType>& a_rhs) {
  volume_m += a_rhs.volume_m;
  return (*this);
}

template <class ScalarType>
inline VolumeBase<ScalarType>& VolumeBase<ScalarType>::operator*=(
    const ScalarType a_rhs) {
  volume_m *= a_rhs;
  return (*this);
}

template <class ScalarType>
inline VolumeBase<ScalarType>& VolumeBase<ScalarType>::operator/=(
    const ScalarType a_rhs) {
  volume_m /= a_rhs;
  return (*this);
}

template <class ScalarType>
inline VolumeBase<ScalarType>::operator ScalarType() const {
  return volume_m;
}

template <class ScalarType>
inline VolumeBase<ScalarType>& VolumeBase<ScalarType>::operator=(
    const ScalarType a_value) {
  volume_m = a_value;
  return (*this);
}

template <class ScalarType>
inline std::ostream& operator<<(std::ostream& out,
                                const VolumeBase<ScalarType>& a_volume) {
  out << static_cast<ScalarType>(a_volume);
  return out;
}

}  // namespace IRL

#endif  // IRL_MOMENTS_VOLUME_TPP_
