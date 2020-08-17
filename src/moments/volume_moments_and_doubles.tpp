// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_VOLUME_MOMENTS_AND_DOUBLES_TPP_
#define SRC_MOMENTS_VOLUME_MOMENTS_AND_DOUBLES_TPP_

#include <cstring>

namespace IRL {

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength>::VolumeMomentsAndDoubles(void)
    : volume_moments_m(VolumeMoments::fromScalarConstant(0.0)) {
  data_moments_m.fill(0.0);
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength>::VolumeMomentsAndDoubles(
    const double a_volume, const Pt& a_centroid,
    const ArrayType& a_initial_data)
    : volume_moments_m(a_volume, a_centroid) {
  data_moments_m = a_initial_data;
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength>
VolumeMomentsAndDoubles<kArrayLength>::fromRawDoublePointer(
    const double* a_list, const double* a_data_list) {
  return VolumeMomentsAndDoubles<kArrayLength>(a_list, a_data_list);
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength>
VolumeMomentsAndDoubles<kArrayLength>::fromScalarConstant(
    const double a_value) {
  return VolumeMomentsAndDoubles<kArrayLength>(a_value, a_value);
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength>
VolumeMomentsAndDoubles<kArrayLength>::fromScalarConstant(
    const double a_value, const double a_value_for_data) {
  return VolumeMomentsAndDoubles<kArrayLength>(a_value, a_value_for_data);
}

template <UnsignedIndex_t kArrayLength>
template <class GeometryType>
inline VolumeMomentsAndDoubles<kArrayLength>
VolumeMomentsAndDoubles<kArrayLength>::calculateMoments(
    GeometryType* a_geometry) {
  return a_geometry->calculateVolumeMomentsAndDoubles();
}

template <UnsignedIndex_t kArrayLength>
inline Volume& VolumeMomentsAndDoubles<kArrayLength>::volume(void) {
  return volume_moments_m.volume();
}

template <UnsignedIndex_t kArrayLength>
inline const Volume& VolumeMomentsAndDoubles<kArrayLength>::volume(void) const {
  return volume_moments_m.volume();
}

template <UnsignedIndex_t kArrayLength>
inline Pt& VolumeMomentsAndDoubles<kArrayLength>::centroid(void) {
  return volume_moments_m.centroid();
}

template <UnsignedIndex_t kArrayLength>
inline const Pt& VolumeMomentsAndDoubles<kArrayLength>::centroid(void) const {
  return volume_moments_m.centroid();
}

template <UnsignedIndex_t kArrayLength>
inline typename VolumeMomentsAndDoubles<kArrayLength>::ArrayType&
VolumeMomentsAndDoubles<kArrayLength>::data(void) {
  return data_moments_m;
}

template <UnsignedIndex_t kArrayLength>
inline const typename VolumeMomentsAndDoubles<kArrayLength>::ArrayType&
VolumeMomentsAndDoubles<kArrayLength>::data(void) const {
  return data_moments_m;
}

template <UnsignedIndex_t kArrayLength>
inline void VolumeMomentsAndDoubles<kArrayLength>::normalizeByVolume(void) {
  volume_moments_m.normalizeByVolume();
  if (this->volume() != 0.0) {
    for (UnsignedIndex_t n = 0; n < data_moments_m.size(); ++n) {
      this->data()[n] /= this->volume();
    }
  }
}

template <UnsignedIndex_t kArrayLength>
inline void VolumeMomentsAndDoubles<kArrayLength>::multiplyByVolume(void) {
  volume_moments_m.multiplyByVolume();
  for (UnsignedIndex_t n = 0; n < data_moments_m.size(); ++n) {
    this->data()[n] *= this->volume();
  }
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength>&
VolumeMomentsAndDoubles<kArrayLength>::operator+=(
    const VolumeMomentsAndDoubles<kArrayLength>& a_rhs) {
  volume_moments_m += a_rhs.volume_moments_m;
  for (UnsignedIndex_t n = 0; n < data_moments_m.size(); ++n) {
    this->data()[n] += a_rhs.data()[n];
  }
  return (*this);
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength>&
VolumeMomentsAndDoubles<kArrayLength>::operator*=(const double a_rhs) {
  volume_moments_m *= a_rhs;
  for (UnsignedIndex_t n = 0; n < data_moments_m.size(); ++n) {
    this->data()[n] *= a_rhs;
  }
  return (*this);
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength>&
VolumeMomentsAndDoubles<kArrayLength>::operator/=(const double a_rhs) {
  volume_moments_m /= a_rhs;
  for (UnsignedIndex_t n = 0; n < data_moments_m.size(); ++n) {
    this->data()[n] /= a_rhs;
  }
  return (*this);
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength>&
VolumeMomentsAndDoubles<kArrayLength>::operator=(const double a_value) {
  volume_moments_m = a_value;
  for (UnsignedIndex_t n = 0; n < data_moments_m.size(); ++n) {
    this->data()[n] = a_value;
  }
  return (*this);
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength>::VolumeMomentsAndDoubles(
    const double* a_list, const double* a_data_list)
    : volume_moments_m(VolumeMoments::fromRawDoublePointer(a_list)) {
  std::memcpy(data_moments_m.data(), a_data_list,
              kArrayLength * sizeof(double));
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength>::VolumeMomentsAndDoubles(
    const double a_value, const double a_value_for_data)
    : volume_moments_m(VolumeMoments::fromScalarConstant(a_value)) {
  data_moments_m.fill(a_value_for_data);
}

template <UnsignedIndex_t kArrayLength>
inline std::ostream& operator<<(
    std::ostream& out,
    const VolumeMomentsAndDoubles<kArrayLength>& a_volume_moments) {
  out << static_cast<double>(a_volume_moments.volume()) << " ";
  out << a_volume_moments.centroid();
  out << " Data: ";
  for (UnsignedIndex_t n = 0; n < kArrayLength; ++n) {
    out << a_volume_moments.data()[n] << " ";
  }
  return out;
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength> operator+(
    const VolumeMomentsAndDoubles<kArrayLength>& a_vmad1,
    const VolumeMomentsAndDoubles<kArrayLength>& a_vmad2) {
  VolumeMomentsAndDoubles<kArrayLength> vmad_to_return;
  vmad_to_return.volume() = a_vmad1.volume() + a_vmad2.volume();
  vmad_to_return.centroid() = a_vmad1.centroid() + a_vmad2.centroid();
  for (UnsignedIndex_t n = 0; n < kArrayLength; ++n) {
    vmad_to_return.data()[n] = a_vmad1.data()[n] + a_vmad2.data()[n];
  }
  return vmad_to_return;
}
template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength> operator-(
    const VolumeMomentsAndDoubles<kArrayLength>& a_vmad1,
    const VolumeMomentsAndDoubles<kArrayLength>& a_vmad2) {
  VolumeMomentsAndDoubles<kArrayLength> vmad_to_return;
  vmad_to_return.volume() = a_vmad1.volume() - a_vmad2.volume();
  vmad_to_return.centroid() = a_vmad1.centroid() - a_vmad2.centroid();
  for (UnsignedIndex_t n = 0; n < kArrayLength; ++n) {
    vmad_to_return.data()[n] = a_vmad1.data()[n] - a_vmad2.data()[n];
  }
  return vmad_to_return;
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength> operator*(
    const double a_multiplier,
    const VolumeMomentsAndDoubles<kArrayLength>& a_vmad) {
  VolumeMomentsAndDoubles<kArrayLength> vmad_to_return;
  vmad_to_return.volume() = a_vmad.volume() * a_multiplier;
  vmad_to_return.centroid() = a_vmad.centroid() * a_multiplier;
  for (UnsignedIndex_t n = 0; n < kArrayLength; ++n) {
    vmad_to_return.data()[n] = a_vmad.data()[n] * a_multiplier;
  }
  return vmad_to_return;
}

template <UnsignedIndex_t kArrayLength>
inline VolumeMomentsAndDoubles<kArrayLength> operator*(
    const VolumeMomentsAndDoubles<kArrayLength>& a_vmad,
    const double a_multiplier) {
  return a_multiplier * a_vmad;
}

}  // namespace IRL

#endif  // SRC_MOMENTS_VOLUME_MOMENTS_AND_DOUBLES_TPP_
