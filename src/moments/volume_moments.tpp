// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_VOLUME_MOMENTS_TPP_
#define SRC_MOMENTS_VOLUME_MOMENTS_TPP_

namespace IRL {

inline VolumeMoments::VolumeMoments(void)
    : volume_m{0.0}, centroid_m{0.0, 0.0, 0.0} {}

inline constexpr VolumeMoments::VolumeMoments(const double a_volume,
                                              const Pt& a_centroid)
    : volume_m(a_volume), centroid_m(a_centroid) {}

inline constexpr VolumeMoments VolumeMoments::fromRawDoublePointer(
    const double* a_list) {
  return VolumeMoments(a_list);
}

inline constexpr VolumeMoments VolumeMoments::fromScalarConstant(
    const double a_value) {
  return VolumeMoments(a_value);
}

template <class GeometryType>
inline VolumeMoments VolumeMoments::calculateMoments(GeometryType* a_geometry) {
  return a_geometry->calculateMoments();
}

inline Volume& VolumeMoments::volume(void) { return volume_m; }

inline const Volume& VolumeMoments::volume(void) const { return volume_m; }

inline Pt& VolumeMoments::centroid(void) { return centroid_m; }

inline const Pt& VolumeMoments::centroid(void) const { return centroid_m; }

inline void VolumeMoments::normalizeByVolume(void) {
  if (this->volume() != 0.0) {
    this->centroid() /= this->volume();
  }
}

inline void VolumeMoments::multiplyByVolume(void) {
  this->centroid() *= this->volume();
}

inline VolumeMoments& VolumeMoments::operator+=(const VolumeMoments& a_rhs) {
  this->volume() += a_rhs.volume();
  this->centroid() += a_rhs.centroid();
  return (*this);
}

inline VolumeMoments& VolumeMoments::operator*=(const double a_rhs) {
  this->volume() *= a_rhs;
  this->centroid() *= a_rhs;
  return (*this);
}

inline VolumeMoments& VolumeMoments::operator/=(const double a_rhs) {
  this->volume() /= a_rhs;
  this->centroid() /= a_rhs;
  return (*this);
}

inline VolumeMoments& VolumeMoments::operator=(const double a_value) {
  this->volume() = a_value;
  this->centroid() = a_value;
  return (*this);
}

inline constexpr VolumeMoments::VolumeMoments(const double* a_list)
    : volume_m(a_list[0]), centroid_m{a_list[1], a_list[2], a_list[3]} {}

inline constexpr VolumeMoments::VolumeMoments(const double a_value)
    : volume_m(a_value), centroid_m(Pt::fromScalarConstant(a_value)) {}

inline std::ostream& operator<<(std::ostream& out,
                                const VolumeMoments& a_volume_moments) {
  out << static_cast<double>(a_volume_moments.volume()) << " ";
  out << a_volume_moments.centroid();
  return out;
}

inline VolumeMoments operator+(const VolumeMoments& a_vm1,
                               const VolumeMoments& a_vm2) {
  return VolumeMoments(a_vm1.volume() + a_vm2.volume(),
                       IRL::Pt(a_vm1.centroid() + a_vm2.centroid()));
}
inline VolumeMoments operator-(const VolumeMoments& a_vm1,
                               const VolumeMoments& a_vm2) {
  return VolumeMoments(a_vm1.volume() - a_vm2.volume(),
                       a_vm1.centroid() - a_vm2.centroid());
}

inline VolumeMoments operator*(const double a_multiplier,
                               const VolumeMoments& a_vm) {
  return VolumeMoments(a_multiplier * a_vm.volume(),
                       a_multiplier * a_vm.centroid());
}
inline VolumeMoments operator*(const VolumeMoments& a_vm,
                               const double a_multiplier) {
  return a_multiplier * a_vm;
}

}  // namespace IRL

#endif  // SRC_MOMENTS_VOLUME_MOMENTS_TPP_
