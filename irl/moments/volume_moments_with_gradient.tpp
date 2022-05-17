// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_VOLUME_MOMENTS_WITH_GRADIENT_TPP_
#define IRL_MOMENTS_VOLUME_MOMENTS_WITH_GRADIENT_TPP_

namespace IRL {

template <class GradientType>
inline VolumeMomentsWithGradient<GradientType>::VolumeMomentsWithGradient(
    void) {
  volume_m = 0.0;
  volume_gradient_m = GradientType(0.0);
  centroid_m = PtWithGradient<GradientType>();
}

template <class GradientType>
inline constexpr VolumeMomentsWithGradient<
    GradientType>::VolumeMomentsWithGradient(const double a_value) {
  volume_m = a_value;
  volume_gradient_m = GradientType(0.0);
  centroid_m = PtWithGradient<GradientType>();
}
template <class GradientType>

inline constexpr VolumeMomentsWithGradient<GradientType>::
    VolumeMomentsWithGradient(const double a_volume,
                              const GradientType& a_volume_gradient,
                              const PtWithGradient<GradientType>& a_centroid) {
  volume_m = a_volume;
  volume_gradient_m = a_volume_gradient;
  centroid_m = a_centroid;
}

template <class GradientType>
inline VolumeMomentsWithGradient<GradientType>
VolumeMomentsWithGradient<GradientType>::fromScalarConstant(
    const double a_value) {
  return VolumeMomentsWithGradient<GradientType>(a_value);
}

template <class GradientType>
inline void VolumeMomentsWithGradient<GradientType>::multiplyByVolume(void) {
  centroid_m *= volume_m;
}

template <class GradientType>
inline void VolumeMomentsWithGradient<GradientType>::normalizeByVolume(void) {
  centroid_m /= safelyEpsilon(volume_m);
}

template <class GradientType>
template <class GeometryType>
inline VolumeMomentsWithGradient<GradientType>
VolumeMomentsWithGradient<GradientType>::calculateMoments(
    GeometryType* a_geometry) {
  return VolumeMomentsWithGradient<GradientType>(a_geometry->calculateVolume());
}
template <class GradientType>
inline Volume& VolumeMomentsWithGradient<GradientType>::volume(void) {
  return volume_m;
}

template <class GradientType>
inline const Volume& VolumeMomentsWithGradient<GradientType>::volume(
    void) const {
  return volume_m;
}
template <class GradientType>
inline PtWithGradient<GradientType>&
VolumeMomentsWithGradient<GradientType>::centroid(void) {
  return centroid_m;
}

template <class GradientType>
inline const PtWithGradient<GradientType>&
VolumeMomentsWithGradient<GradientType>::centroid(void) const {
  return centroid_m;
}
template <class GradientType>
inline GradientType& VolumeMomentsWithGradient<GradientType>::volume_gradient(
    void) {
  return volume_gradient_m;
}
template <class GradientType>
inline const GradientType&
VolumeMomentsWithGradient<GradientType>::volume_gradient(void) const {
  return volume_gradient_m;
}

template <class GradientType>
inline VolumeMomentsWithGradient<GradientType>&
VolumeMomentsWithGradient<GradientType>::operator+=(
    const VolumeMomentsWithGradient<GradientType>& a_rhs) {
  volume_m += a_rhs.volume_m;
  volume_gradient_m += a_rhs.volume_gradient_m;
  centroid_m += a_rhs.centroid_m;
  return (*this);
}

template <class GradientType>
inline VolumeMomentsWithGradient<GradientType>&
VolumeMomentsWithGradient<GradientType>::operator*=(const double a_rhs) {
  volume_m *= a_rhs;
  volume_gradient_m *= a_rhs;
  centroid_m *= a_rhs;
  return (*this);
}

template <class GradientType>
inline VolumeMomentsWithGradient<GradientType>&
VolumeMomentsWithGradient<GradientType>::operator/=(const double a_rhs) {
  volume_m /= a_rhs;
  volume_gradient_m /= a_rhs;
  centroid_m /= a_rhs;
  return (*this);
}

template <class GradientType>
inline VolumeMomentsWithGradient<GradientType> operator+(
    const VolumeMomentsWithGradient<GradientType>& a_vm1,
    const VolumeMomentsWithGradient<GradientType>& a_vm2) {
  return VolumeMomentsWithGradient<GradientType>(
      a_vm1.volume() + a_vm2.volume(),
      a_vm1.volume_gradient() + a_vm2.volume_gradient(),
      a_vm1.centroid() + a_vm2.centroid());
}

template <class GradientType>
inline VolumeMomentsWithGradient<GradientType> operator-(
    const VolumeMomentsWithGradient<GradientType>& a_vm1,
    const VolumeMomentsWithGradient<GradientType>& a_vm2) {
  return VolumeMomentsWithGradient<GradientType>(
      a_vm1.volume() - a_vm2.volume(),
      a_vm1.volume_gradient() - a_vm2.volume_gradient(),
      a_vm1.centroid() - a_vm2.centroid());
}

template <class GradientType>
inline VolumeMomentsWithGradient<GradientType> operator*(
    const double a_multiplier,
    const VolumeMomentsWithGradient<GradientType>& a_vm) {
  return VolumeMomentsWithGradient<GradientType>(
      a_multiplier * a_vm.volume(), a_multiplier * a_vm.volume_gradient(),
      a_multiplier * a_vm.centroid());
}

template <class GradientType>
inline VolumeMomentsWithGradient<GradientType> operator*(
    const VolumeMomentsWithGradient<GradientType>& a_vm,
    const double a_multiplier) {
  return VolumeMomentsWithGradient<GradientType>(
      a_multiplier * a_vm.volume(), a_multiplier * a_vm.volume_gradient(),
      a_multiplier * a_vm.centroid());
}

}  // namespace IRL

#endif  // IRL_MOMENTS_VOLUME_MOMENTS_WITH_GRADIENT_TPP_
