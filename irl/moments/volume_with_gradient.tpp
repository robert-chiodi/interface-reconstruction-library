// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_VOLUME_WITH_GRADIENT_TPP_
#define IRL_MOMENTS_VOLUME_WITH_GRADIENT_TPP_

namespace IRL {

template <class GradientType>
inline VolumeWithGradient<GradientType>::VolumeWithGradient(void)
    : volume_m{0.0}, gradient_m{0.0} {}

template <class GradientType>
inline constexpr VolumeWithGradient<GradientType>::VolumeWithGradient(
    const double a_value)
    : volume_m{a_value}, gradient_m{0.0} {}

template <class GradientType>
inline VolumeWithGradient<GradientType>
VolumeWithGradient<GradientType>::fromScalarConstant(const double a_value) {
  return VolumeWithGradient<GradientType>(a_value);
}

template <class GradientType>
inline void VolumeWithGradient<GradientType>::multiplyByVolume(void) {}

template <class GradientType>
inline void VolumeWithGradient<GradientType>::normalizeByVolume(void) {}

template <class GradientType>
template <class GeometryType>
inline VolumeWithGradient<GradientType>
VolumeWithGradient<GradientType>::calculateMoments(GeometryType* a_geometry) {
  return VolumeWithGradient<GradientType>();
  //   return a_geometry->calculateVolumeAndGradient();
}
template <class GradientType>
inline Volume& VolumeWithGradient<GradientType>::volume(void) {
  return volume_m;
}

template <class GradientType>
inline const Volume& VolumeWithGradient<GradientType>::volume(void) const {
  return volume_m;
}

template <class GradientType>
inline VolumeWithGradient<GradientType>&
VolumeWithGradient<GradientType>::operator+=(
    const VolumeWithGradient<GradientType>& a_rhs) {
  volume_m += a_rhs.volume_m;
  gradient_m += a_rhs.gradient_m;
  return (*this);
}

template <class GradientType>
inline VolumeWithGradient<GradientType>&
VolumeWithGradient<GradientType>::operator*=(const double a_rhs) {
  volume_m *= a_rhs;
  gradient_m *= a_rhs;
  return (*this);
}

template <class GradientType>
inline VolumeWithGradient<GradientType>&
VolumeWithGradient<GradientType>::operator/=(const double a_rhs) {
  volume_m /= a_rhs;
  gradient_m /= a_rhs;
  return (*this);
}

template <class GradientType>
inline VolumeWithGradient<GradientType>::operator double() const {
  return volume_m;
}

template <class GradientType>
inline VolumeWithGradient<GradientType>&
VolumeWithGradient<GradientType>::operator=(const double a_value) {
  volume_m = a_value;
  return (*this);
}

template <class GradientType>
void VolumeWithGradient<GradientType>::isVolumeWithGradient(void) {}

}  // namespace IRL

#endif  // IRL_MOMENTS_VOLUME_WITH_GRADIENT_TPP_
