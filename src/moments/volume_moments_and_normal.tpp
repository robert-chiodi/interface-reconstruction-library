// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_VOLUME_MOMENTS_AND_NORMAL_TPP_
#define SRC_MOMENTS_VOLUME_MOMENTS_AND_NORMAL_TPP_

namespace IRL {

inline VolumeMomentsAndNormal::VolumeMomentsAndNormal(
    const VolumeMoments& a_volume_moments, const Normal& a_normal)
    : volume_moments_m(a_volume_moments), normal_m(a_normal) {}

template <class GeometryType>
inline VolumeMomentsAndNormal VolumeMomentsAndNormal::calculateMoments(
    GeometryType* a_geometry) {
  return a_geometry->calculateVolumeMomentsAndNormal();
}

inline VolumeMomentsAndNormal VolumeMomentsAndNormal::fromScalarConstant(
    const double a_constant) {
  return VolumeMomentsAndNormal(a_constant);
}

inline VolumeMoments& VolumeMomentsAndNormal::volumeMoments(void) {
  return volume_moments_m;
}

inline const VolumeMoments& VolumeMomentsAndNormal::volumeMoments(void) const {
  return volume_moments_m;
}

inline Normal& VolumeMomentsAndNormal::normal(void) { return normal_m; }

inline const Normal& VolumeMomentsAndNormal::normal(void) const {
  return normal_m;
}

inline void VolumeMomentsAndNormal::normalizeByVolume(void) {
  this->volumeMoments().normalizeByVolume();
  if (this->volumeMoments().volume() != 0.0) {
    this->normal() /= this->volumeMoments().volume();
  }
}

inline void VolumeMomentsAndNormal::multiplyByVolume(void) {
  this->volumeMoments().multiplyByVolume();
  this->normal() *= this->volumeMoments().volume();
}

inline VolumeMomentsAndNormal& VolumeMomentsAndNormal::operator+=(
    const VolumeMomentsAndNormal& a_rhs) {
  this->volumeMoments() += a_rhs.volumeMoments();
  this->normal() += a_rhs.normal();
  return (*this);
}

inline VolumeMomentsAndNormal& VolumeMomentsAndNormal::operator*=(
    const double& a_rhs) {
  this->volumeMoments() *= a_rhs;
  this->normal() *= a_rhs;
  return (*this);
}

inline VolumeMomentsAndNormal& VolumeMomentsAndNormal::operator/=(
    const double& a_rhs) {
  this->volumeMoments() /= a_rhs;
  this->normal() /= a_rhs;
  return (*this);
}

inline VolumeMomentsAndNormal& VolumeMomentsAndNormal::operator=(
    const double a_value) {
  this->volumeMoments() = a_value;
  this->normal() = a_value;
  return (*this);
}

inline VolumeMomentsAndNormal::VolumeMomentsAndNormal(const double a_constant)
    : volume_moments_m(VolumeMoments::fromScalarConstant(a_constant)),
      normal_m(Normal::fromScalarConstant(a_constant)) {}

inline std::ostream& operator<<(
    std::ostream& out,
    const VolumeMomentsAndNormal& a_volume_moments_and_normal) {
  out << '\n';
  out << "VolumeMoments : " << a_volume_moments_and_normal.volumeMoments()
      << '\n';
  out << "Normal vector : " << a_volume_moments_and_normal.normal() << '\n';
  return out;
}
}  // namespace IRL

#endif  // SRC_MOMENTS_VOLUME_MOMENTS_AND_NORMAL_TPP_
