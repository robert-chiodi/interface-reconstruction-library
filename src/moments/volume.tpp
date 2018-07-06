// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_VOLUME_TPP_
#define SRC_MOMENTS_VOLUME_TPP_


namespace IRL {

inline Volume::Volume(void) : volume_m{0.0} {}

inline constexpr Volume::Volume(const double a_value) : volume_m{a_value} {}

inline Volume Volume::fromScalarConstant(const double a_value) {
  return Volume(a_value);
}

template <class GeometryType>
inline Volume Volume::calculateMoments(GeometryType* a_geometry) {
  return a_geometry->calculateVolume();
}

inline void Volume::multiplyByVolume(void) {}

inline void Volume::normalizeByVolume(void) {}

inline Volume& Volume::operator+=(const Volume& a_rhs) {
  volume_m += a_rhs.volume_m;
  return (*this);
}

inline Volume& Volume::operator*=(const double a_rhs) {
  volume_m *= a_rhs;
  return (*this);
}

inline Volume& Volume::operator/=(const double a_rhs) {
  volume_m /= a_rhs;
  return (*this);
}

inline Volume::operator double() const { return volume_m; }

inline Volume& Volume::operator=(const double a_value) {
  volume_m = a_value;
  return (*this);
}

inline std::ostream& operator<<(std::ostream& out, const Volume& a_volume) {
  out << static_cast<double>(a_volume);
  return out;
}

}  // namespace IRL

#endif  // SRC_MOMENTS_VOLUME_TPP_
