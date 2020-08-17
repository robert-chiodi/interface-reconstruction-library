// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_SEPARATED_VOLUME_MOMENTS_TPP_
#define SRC_MOMENTS_SEPARATED_VOLUME_MOMENTS_TPP_

namespace IRL {

template <class Derived, class MomentsType>
Derived& SeparatedMomentsCommon<Derived, MomentsType>::getDerived(void) {
  return static_cast<Derived&>(*this);
}
template <class Derived, class MomentsType>
const Derived& SeparatedMomentsCommon<Derived, MomentsType>::getDerived(
    void) const {
  return static_cast<const Derived&>(*this);
}

template <class Derived, class MomentsType>
inline constexpr SeparatedMomentsCommon<
    Derived, MomentsType>::SeparatedMomentsCommon(const MomentsType& a_moment_0,
                                                  const MomentsType& a_moment_1)
    : volume_moments_m{a_moment_0, a_moment_1} {}

template <class Derived, class MomentsType>
inline Derived
SeparatedMomentsCommon<Derived, MomentsType>::fromRawDoublePointer(
    const double* a_list) {
  return Derived(a_list);
}

template <class Derived, class MomentsType>
inline Derived SeparatedMomentsCommon<Derived, MomentsType>::fromScalarConstant(
    const double a_value) {
  return Derived(a_value);
}

template <class Derived, class MomentsType>
template <class GeometryType>
inline Derived
SeparatedMomentsCommon<Derived, MomentsType>::fillWithComplementMoments(
    const MomentsType& a_known_moments,
    const GeometryType& a_encompassing_geometry, const bool a_flipped) {
  return Derived::fillWithComplementMoments(a_known_moments,
                                            a_encompassing_geometry, a_flipped);
}

template <class Derived, class MomentsType>
inline Derived
SeparatedMomentsCommon<Derived, MomentsType>::fillWithComplementMoments(
    const MomentsType& a_known_moments,
    const MomentsType& a_encompassing_geometry_volume_moments,
    const bool a_flipped) {
  return Derived::fillWithComplementMoments(
      a_known_moments, a_encompassing_geometry_volume_moments, a_flipped);
}

template <class Derived, class MomentsType>
inline constexpr UnsignedIndex_t
SeparatedMomentsCommon<Derived, MomentsType>::getNumberOfPhases(void) {
  return kMaxNumberOfPhases;
}

template <class Derived, class MomentsType>
inline MomentsType& SeparatedMomentsCommon<Derived, MomentsType>::operator[](
    const UnsignedIndex_t a_moment_index) {
  assert(a_moment_index < this->getNumberOfPhases());
  return volume_moments_m[a_moment_index];
}

template <class Derived, class MomentsType>
inline const MomentsType& SeparatedMomentsCommon<Derived, MomentsType>::
operator[](const UnsignedIndex_t a_moment_index) const {
  assert(a_moment_index < this->getNumberOfPhases());
  return volume_moments_m[a_moment_index];
}

template <class Derived, class MomentsType>
inline Derived& SeparatedMomentsCommon<Derived, MomentsType>::operator+=(
    const Derived& a_rhs) {
  assert(this->getNumberOfPhases() == a_rhs.getNumberOfPhases());
  for (UnsignedIndex_t n = 0; n < this->getNumberOfPhases(); ++n) {
    volume_moments_m[n] += a_rhs[n];
  }
  return this->getDerived();
}

template <class Derived, class MomentsType>
inline Derived& SeparatedMomentsCommon<Derived, MomentsType>::operator*=(
    const double a_rhs) {
  for (auto& elem : volume_moments_m) {
    elem *= a_rhs;
  }
  return this->getDerived();
}

template <class Derived, class MomentsType>
inline Derived& SeparatedMomentsCommon<Derived, MomentsType>::operator=(
    const double a_value) {
  for (auto& elem : volume_moments_m) {
    elem = a_value;
  }
  return this->getDerived();
}

template <class Derived, class MomentsType>
inline void SeparatedMomentsCommon<Derived, MomentsType>::normalizeByVolume(
    void) {
  for (auto& elem : volume_moments_m) {
    elem.normalizeByVolume();
  }
}

template <class Derived, class MomentsType>
inline void SeparatedMomentsCommon<Derived, MomentsType>::multiplyByVolume(
    void) {
  for (auto& elem : volume_moments_m) {
    elem.multiplyByVolume();
  }
}

template <class Derived, class MomentsType>
inline void SeparatedMomentsCommon<Derived, MomentsType>::swap(
    const UnsignedIndex_t a_index_0, const UnsignedIndex_t a_index_1) {
  assert(a_index_0 < this->getNumberOfPhases());
  assert(a_index_1 < this->getNumberOfPhases());
  std::swap((*this)[a_index_0], (*this)[a_index_1]);
}

template <class Derived, class MomentsType>
inline typename SeparatedMomentsCommon<Derived, MomentsType>::iterator
SeparatedMomentsCommon<Derived, MomentsType>::begin(void) noexcept {
  return volume_moments_m.begin();
}
template <class Derived, class MomentsType>
inline typename SeparatedMomentsCommon<Derived, MomentsType>::const_iterator
SeparatedMomentsCommon<Derived, MomentsType>::begin(void) const noexcept {
  return this->cbegin();
}
template <class Derived, class MomentsType>
inline typename SeparatedMomentsCommon<Derived, MomentsType>::const_iterator
SeparatedMomentsCommon<Derived, MomentsType>::cbegin(void) const noexcept {
  return volume_moments_m.cbegin();
}
template <class Derived, class MomentsType>
inline typename SeparatedMomentsCommon<Derived, MomentsType>::iterator
SeparatedMomentsCommon<Derived, MomentsType>::end(void) noexcept {
  return volume_moments_m.end();
}
template <class Derived, class MomentsType>
inline typename SeparatedMomentsCommon<Derived, MomentsType>::const_iterator
SeparatedMomentsCommon<Derived, MomentsType>::end(void) const noexcept {
  return this->cend();
}
template <class Derived, class MomentsType>
inline typename SeparatedMomentsCommon<Derived, MomentsType>::const_iterator
SeparatedMomentsCommon<Derived, MomentsType>::cend(void) const noexcept {
  return volume_moments_m.cend();
}

template <class Derived, class MomentsType>
inline SeparatedMomentsCommon<Derived, MomentsType>::SeparatedMomentsCommon(
    const double a_value)
    : volume_moments_m{MomentsType::fromScalarConstant(a_value),
                       MomentsType::fromScalarConstant(a_value)} {}

template <class GeometryType>
inline SeparatedMoments<VolumeMoments>
SeparatedMoments<VolumeMoments>::fillWithComplementMoments(
    const VolumeMoments& a_known_moments,
    const GeometryType& a_encompassing_geometry, const bool a_flipped) {
  // Here, we are expecting a_known_moments.centroid = vol*centroid, and
  // NOT the centroid. Will return correct centroids for both phases
  MomentsType encompassing_moments = a_encompassing_geometry.calculateMoments();

  // Encompassing volume should be same or greater than the known phase's
  // volume. Doing for round-off
  MomentsType unknown_moments;
//  if (std::fabs(encompassing_moments.volume()) <
//      std::fabs(a_known_moments.volume())) {
//    unknown_moments = MomentsType::fromScalarConstant(0.0);
//  } else {
    unknown_moments = encompassing_moments - a_known_moments;
//  }
  return a_flipped ? SelfType(unknown_moments, a_known_moments)
                   : SelfType(a_known_moments, unknown_moments);
}

inline SeparatedMoments<VolumeMoments>
SeparatedMoments<VolumeMoments>::fillWithComplementMoments(
    const VolumeMoments& a_known_moments,
    const VolumeMoments& a_encompassing_geometry_volume_moments,
    const bool a_flipped) {
  // Here, we are expecting a_encompassing_geometry_volume_moments.centroid =
  // vol*centroid, and NOT the centroid. Will return correct centroids for
  // both phases

  // Encompassing volume should be same or greater than the known phase's
  // volume.  Doing for round-off
  MomentsType unknown_moments;
//  if (std::fabs(a_encompassing_geometry_volume_moments.volume()) <
//      std::fabs(a_known_moments.volume())) {
//    unknown_moments = MomentsType::fromScalarConstant(0.0);
//  } else {
    unknown_moments = a_encompassing_geometry_volume_moments - a_known_moments;
//  }

  return a_flipped ? SelfType(unknown_moments, a_known_moments)
                   : SelfType(a_known_moments, unknown_moments);
}

inline SeparatedMoments<VolumeMoments>::SeparatedMoments(const double* a_list)
    : SeparatedMomentsCommon<SeparatedMoments<VolumeMoments>,MomentsType>(MomentsType::fromRawDoublePointer(&a_list[0]),
                       MomentsType::fromRawDoublePointer(&a_list[4])) {}


template <class GeometryType>
inline SeparatedMoments<Volume>
SeparatedMoments<Volume>::fillWithComplementMoments(
    const MomentsType& a_known_moments,
    const GeometryType& a_encompassing_geometry, const bool a_flipped) {

  MomentsType encompassing_moments = a_encompassing_geometry.calculateVolume();

  // Encompassing volume should be same or greater than the known phase's
  // volume. Doing for round-off
  MomentsType unknown_moments;
//  if (std::fabs(encompassing_moments) <
//      std::fabs(a_known_moments)) {
//    unknown_moments = MomentsType::fromScalarConstant(0.0);
//  } else {
    unknown_moments = encompassing_moments - a_known_moments;
//  }
  return a_flipped ? SelfType(unknown_moments, a_known_moments)
                   : SelfType(a_known_moments, unknown_moments);
}

inline SeparatedMoments<Volume>
SeparatedMoments<Volume>::fillWithComplementMoments(
    const MomentsType& a_known_moments,
    const MomentsType& a_encompassing_geometry_volume_moments,
    const bool a_flipped) {

  // Encompassing volume should be same or greater than the known phase's
  // volume.  Doing for round-off
  MomentsType unknown_moments;
//  if (std::fabs(a_encompassing_geometry_volume_moments) <
//      std::fabs(a_known_moments)) {
//    unknown_moments = MomentsType::fromScalarConstant(0.0);
//  } else {
    unknown_moments = a_encompassing_geometry_volume_moments - a_known_moments;
//  }

  return a_flipped ? SelfType(unknown_moments, a_known_moments)
                   : SelfType(a_known_moments, unknown_moments);
}

inline SeparatedMoments<Volume>::SeparatedMoments(const double* a_list)
    : SeparatedMomentsCommon<SeparatedMoments<Volume>,MomentsType>(MomentsType(a_list[0]),
                                                                   MomentsType(a_list[1])) {}



template <UnsignedIndex_t kArrayLength>
template <class GeometryType>
inline SeparatedMoments<VolumeMomentsAndDoubles<kArrayLength>>
SeparatedMoments<VolumeMomentsAndDoubles<kArrayLength>>::
    fillWithComplementMoments(
        const VolumeMomentsAndDoubles<kArrayLength>& a_known_moments,
        const GeometryType& a_encompassing_geometry, const bool a_flipped) {
  // Here, we are expecting a_known_moments.centroid = vol*centroid, and
  // NOT the centroid. Will return correct centroids for both phases
  MomentsType encompassing_moments =
      MomentsType::calculateMoments(&a_encompassing_geometry);

  // Encompassing volume should be same or greater than the known phase's
  // volume
  MomentsType unknown_moments;
  if (std::fabs(encompassing_moments.volume()) <
      std::fabs(a_known_moments.volume())) {
    unknown_moments = MomentsType::fromScalarConstant(0.0);
  } else {
    unknown_moments = encompassing_moments - a_known_moments;
  }

  return a_flipped ? SelfType(unknown_moments, a_known_moments)
                   : SelfType(a_known_moments, unknown_moments);
}

template <UnsignedIndex_t kArrayLength>
inline SeparatedMoments<VolumeMomentsAndDoubles<kArrayLength>>
SeparatedMoments<VolumeMomentsAndDoubles<kArrayLength>>::
    fillWithComplementMoments(
        const VolumeMomentsAndDoubles<kArrayLength>& a_known_moments,
        const VolumeMomentsAndDoubles<kArrayLength>&
            a_encompassing_geometry_volume_moments,
        const bool a_flipped) {
  MomentsType unknown_moments;
  if (std::fabs(a_encompassing_geometry_volume_moments.volume()) <
      std::fabs(a_known_moments.volume())) {
    unknown_moments = MomentsType::fromScalarConstant(0.0);
  } else {
    unknown_moments = a_encompassing_geometry_volume_moments - a_known_moments;
  }

  return a_flipped ? SelfType(unknown_moments, a_known_moments)
                   : SelfType(a_known_moments, unknown_moments);
}

template <class MomentsType>
inline std::ostream& operator<<(
    std::ostream& out,
    const SeparatedMoments<MomentsType>& a_separated_volume_moments) {
  out << '\n';
  for (UnsignedIndex_t phase = 0;
       phase < a_separated_volume_moments.getNumberOfPhases(); ++phase) {
    out << "Moments for Phase " << phase << ": "
        << a_separated_volume_moments[phase] << '\n';
  }
  return out;
}

template <class Derived, class MomentsType>
inline SeparatedMomentsCommon<Derived, MomentsType> operator*(
    const SeparatedMomentsCommon<Derived, MomentsType>& a_svm,
    const double a_multiplier) {
  return a_multiplier * a_svm;
}

template <class Derived, class MomentsType>
inline SeparatedMomentsCommon<Derived, MomentsType> operator*(
    const double a_multiplier,
    const SeparatedMomentsCommon<Derived, MomentsType>& a_svm) {
  return Derived(a_multiplier * a_svm[0], a_multiplier * a_svm[1]);
}

}  // namespace IRL

#endif  // SRC_MOMENTS_SEPARATED_VOLUME_MOMENTS_TPP_
