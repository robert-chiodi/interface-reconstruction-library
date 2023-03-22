// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_VOLUME_MOMENTS_TPP_
#define IRL_MOMENTS_VOLUME_MOMENTS_TPP_

namespace IRL {

template <class ScalarType>
inline VolumeMomentsBase<ScalarType>::VolumeMomentsBase(void)
    : volume_m{static_cast<ScalarType>(0)},
      centroid_m{static_cast<ScalarType>(0), static_cast<ScalarType>(0),
                 static_cast<ScalarType>(0)} {}

template <class ScalarType>
inline constexpr VolumeMomentsBase<ScalarType>::VolumeMomentsBase(
    const ScalarType a_volume, const PtBase<ScalarType>& a_centroid)
    : volume_m(a_volume), centroid_m(a_centroid) {}

template <class ScalarType>
inline constexpr VolumeMomentsBase<ScalarType>::VolumeMomentsBase(
    const VolumeMomentsBase<double>& a_moments)
    : volume_m{static_cast<ScalarType>(a_moments.volume())},
      centroid_m{static_cast<ScalarType>(a_moments.centroid()[0]),
                 static_cast<ScalarType>(a_moments.centroid()[1]),
                 static_cast<ScalarType>(a_moments.centroid()[2])} {}

template <class ScalarType>
inline constexpr VolumeMomentsBase<ScalarType>::VolumeMomentsBase(
    const VolumeMomentsBase<Quad_t>& a_moments)
    : volume_m{static_cast<ScalarType>(a_moments.volume())},
      centroid_m{static_cast<ScalarType>(a_moments.centroid()[0]),
                 static_cast<ScalarType>(a_moments.centroid()[1]),
                 static_cast<ScalarType>(a_moments.centroid()[2])} {}

template <class ScalarType>
inline constexpr VolumeMomentsBase<ScalarType>
VolumeMomentsBase<ScalarType>::fromRawDoublePointer(const ScalarType* a_list) {
  return VolumeMomentsBase<ScalarType>(a_list);
}

template <class ScalarType>
inline constexpr VolumeMomentsBase<ScalarType>
VolumeMomentsBase<ScalarType>::fromScalarConstant(const ScalarType a_value) {
  return VolumeMomentsBase<ScalarType>(static_cast<ScalarType>(a_value));
}

template <class ScalarType>
template <class GeometryType>
inline VolumeMomentsBase<ScalarType>
VolumeMomentsBase<ScalarType>::calculateMoments(GeometryType* a_geometry) {
  return VolumeMomentsBase<ScalarType>(a_geometry->calculateMoments());
}

template <class ScalarType>
inline VolumeBase<ScalarType>& VolumeMomentsBase<ScalarType>::volume(void) {
  return volume_m;
}

template <class ScalarType>
inline const VolumeBase<ScalarType>& VolumeMomentsBase<ScalarType>::volume(
    void) const {
  return volume_m;
}

template <class ScalarType>
inline PtBase<ScalarType>& VolumeMomentsBase<ScalarType>::centroid(void) {
  return centroid_m;
}

template <class ScalarType>
inline const PtBase<ScalarType>& VolumeMomentsBase<ScalarType>::centroid(
    void) const {
  return centroid_m;
}

template <class ScalarType>
inline void VolumeMomentsBase<ScalarType>::normalizeByVolume(void) {
  if (this->volume() != static_cast<ScalarType>(0)) {
    this->centroid() /= this->volume();
  }
}

template <class ScalarType>
inline void VolumeMomentsBase<ScalarType>::multiplyByVolume(void) {
  this->centroid() *= this->volume();
}

template <class ScalarType>
inline VolumeMomentsBase<ScalarType>& VolumeMomentsBase<ScalarType>::operator+=(
    const VolumeMomentsBase<ScalarType>& a_rhs) {
  this->volume() += a_rhs.volume();
  this->centroid() += a_rhs.centroid();
  return (*this);
}

template <class ScalarType>
inline VolumeMomentsBase<ScalarType>& VolumeMomentsBase<ScalarType>::operator*=(
    const ScalarType a_rhs) {
  this->volume() *= a_rhs;
  this->centroid() *= a_rhs;
  return (*this);
}

template <class ScalarType>
inline VolumeMomentsBase<ScalarType>& VolumeMomentsBase<ScalarType>::operator/=(
    const ScalarType a_rhs) {
  this->volume() /= a_rhs;
  this->centroid() /= a_rhs;
  return (*this);
}

template <class ScalarType>
inline VolumeMomentsBase<ScalarType>& VolumeMomentsBase<ScalarType>::operator=(
    const ScalarType a_value) {
  this->volume() = a_value;
  this->centroid() = a_value;
  return (*this);
}

template <class ScalarType>
inline VolumeMomentsBase<ScalarType> VolumeMomentsBase<ScalarType>::operator-(
    void) const {
  return VolumeMomentsBase<ScalarType>(-volume_m, -centroid_m);
}

template <class ScalarType>
inline constexpr VolumeMomentsBase<ScalarType>::VolumeMomentsBase(
    const ScalarType* a_list)
    : volume_m(a_list[0]), centroid_m{a_list[1], a_list[2], a_list[3]} {}

template <class ScalarType>
inline constexpr VolumeMomentsBase<ScalarType>::VolumeMomentsBase(
    const ScalarType a_value)
    : volume_m(a_value),
      centroid_m(PtBase<ScalarType>::fromScalarConstant(a_value)) {}

template <class ScalarType>
inline std::ostream& operator<<(
    std::ostream& out, const VolumeMomentsBase<ScalarType>& a_volume_moments) {
  out << static_cast<ScalarType>(a_volume_moments.volume()) << " ";
  out << a_volume_moments.centroid();
  return out;
}

template <class ScalarType>
inline VolumeMomentsBase<ScalarType> operator+(
    const VolumeMomentsBase<ScalarType>& a_vm1,
    const VolumeMomentsBase<ScalarType>& a_vm2) {
  return VolumeMomentsBase<ScalarType>(
      a_vm1.volume() + a_vm2.volume(),
      IRL::PtBase<ScalarType>(a_vm1.centroid() + a_vm2.centroid()));
}

template <class ScalarType>
inline VolumeMomentsBase<ScalarType> operator-(
    const VolumeMomentsBase<ScalarType>& a_vm1,
    const VolumeMomentsBase<ScalarType>& a_vm2) {
  return VolumeMomentsBase<ScalarType>(a_vm1.volume() - a_vm2.volume(),
                                       a_vm1.centroid() - a_vm2.centroid());
}

template <class ScalarType>
inline VolumeMomentsBase<ScalarType> operator*(
    const ScalarType a_multiplier, const VolumeMomentsBase<ScalarType>& a_vm) {
  return VolumeMomentsBase<ScalarType>(a_multiplier * a_vm.volume(),
                                       a_multiplier * a_vm.centroid());
}

template <class ScalarType>
inline VolumeMomentsBase<ScalarType> operator*(
    const VolumeMomentsBase<ScalarType>& a_vm, const ScalarType a_multiplier) {
  return a_multiplier * a_vm;
}

}  // namespace IRL

#endif  // IRL_MOMENTS_VOLUME_MOMENTS_TPP_
