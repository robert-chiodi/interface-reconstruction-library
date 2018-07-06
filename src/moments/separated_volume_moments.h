// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_SEPARATED_VOLUME_MOMENTS_H_
#define SRC_MOMENTS_SEPARATED_VOLUME_MOMENTS_H_

#include <utility>

#include "src/moments/volume_moments.h"
#include "src/moments/volume_moments_and_doubles.h"
#include "src/parameters/defined_types.h"

namespace IRL {
/// \brief Storage for multiple volume moments.
template <class Derived, class MomentsType>
class SeparatedMomentsCommon {
  static constexpr UnsignedIndex_t kMaxNumberOfPhases = 2;
  using iterator = typename std::array<MomentsType, 2>::iterator;
  using const_iterator = typename std::array<MomentsType, 2>::const_iterator;

  Derived& getDerived(void);
  const Derived& getDerived(void) const;

 public:
  using moments_type = MomentsType;

  /// \brief Default constructor.
  SeparatedMomentsCommon(void) = default;

  /// \brief Constructor that sets moments for two phases.
  constexpr SeparatedMomentsCommon(const MomentsType& a_moment_0,
                                   const MomentsType& a_moment_1);

  static Derived fromRawDoublePointer(const double* a_list);

  static Derived fromScalarConstant(const double a_value);

  /// \brief Use knowledge of the initial geometry and moments of the cut phase
  /// to calculate the other phase's moments, return VOLUME WEIGHTED
  /// SeparatedMomentsCommon.
  ///
  /// This function uses knowledge of the initial geomtry (such as
  /// a rectangular cuboid or tet) and the moments known from
  /// the computational cutting to obtain the moments for the other phase.
  ///
  /// NOTICE: The moments in `a_known_moments` are actually volume and
  /// volume*centroid.
  ///
  /// ASIDE: These moments are first assumed to be the liquid phase in
  /// the `SeparatedMomentsCommon` object that is returned. The bool
  /// a_flipped is then use to store in SeparatedMomentsCommon as
  /// {a_known_moments, unknown_moments} (if true) or
  /// {unknown_moments, a_known_moments} if false.
  ///
  /// Template Requirements for `GeomtryType`:
  /// - A `volume(void)` method that returns the volume of the
  /// encompassing geometry
  /// - A `centroid(void)` method that returns the centroid
  ///  (NOT volume*centroid) of the encompassing geometry
  ///
  /// \param[in] a_known_moments GeomtricMoments that are known.
  /// \param[in] a_encompassing_geometry The geomtry from which
  /// the known moments are a subset.
  /// \param[in] a_flipped Boolean expressing whether the supplied
  /// `a_known_moments` is for SeparatedMomentsCommon[0] or [1].
  template <class GeometryType>
  inline static Derived fillWithComplementMoments(
      const MomentsType& a_known_moments,
      const GeometryType& a_encompassing_geometry, const bool a_flipped);

  inline static Derived fillWithComplementMoments(
      const MomentsType& a_known_moments,
      const MomentsType& a_encompassing_geometry_volume_moments,
      const bool a_flipped);

  static constexpr UnsignedIndex_t getNumberOfPhases(void);

  MomentsType& operator[](const UnsignedIndex_t a_moment_index);

  const MomentsType& operator[](const UnsignedIndex_t a_moment_index) const;

  /// \brief Overload += operator to adjust liquid and gas moments.
  Derived& operator+=(const Derived& a_rhs);

  /// \brief Overload *= operator to for const doubles
  Derived& operator*=(const double a_rhs);

  /// \brief Overload assignment to assign constant value to moments.
  Derived& operator=(const double a_value);

  /// \brief Normalize both centroids by their respective phase volumes.
  void normalizeByVolume(void);

  /// \brief Multiply both centroids by their respective phase volumes.
  void multiplyByVolume(void);

  /// \brief Swap two stores moments
  void swap(const UnsignedIndex_t a_index_0, const UnsignedIndex_t a_index_1);

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief  Default destructor.
  ~SeparatedMomentsCommon(void) = default;

 private:

  /// \brief Construct that initializes volume/centroid as a value.
  explicit SeparatedMomentsCommon(const double a_value);

  std::array<MomentsType, kMaxNumberOfPhases>
      volume_moments_m;  ///< \brief VolumeMoments
};

template <class MomentsType>
class SeparatedMoments;

template <>
class SeparatedMoments<VolumeMoments>
    : public SeparatedMomentsCommon<SeparatedMoments<VolumeMoments>,
                                    VolumeMoments> {
  using SelfType = SeparatedMoments<VolumeMoments>;
  using MomentsType = VolumeMoments;

  friend SeparatedMomentsCommon<SelfType,MomentsType>;

 public:
  using SeparatedMomentsCommon<SelfType, MomentsType>::SeparatedMomentsCommon;

  SeparatedMoments(void) = default;

  template <class GeometryType>
  static SeparatedMoments fillWithComplementMoments(
      const MomentsType& a_known_moments,
      const GeometryType& a_encompassing_geometry, const bool a_flipped);

  static SeparatedMoments fillWithComplementMoments(
      const MomentsType& a_known_moments,
      const MomentsType& a_encompassing_geometry_volume_moments,
      const bool a_flipped);

	/// \brief Set volume moments from a list of 8 doubles.
	///
	/// Construct the volume moments from a list of 8 doubles.
	/// Expected order is InternalVolumeMoments(Volume,Centroidx_, Centroid_y,
	/// Centroid_z), then ExternalVolumeMoments(Volume,Centroidx_, Centroid_y,
	/// Centroid_z).
	explicit SeparatedMoments(const double* a_list);


};

template <>
class SeparatedMoments<Volume>
    : public SeparatedMomentsCommon<SeparatedMoments<Volume>,
                                    Volume> {
  using SelfType = SeparatedMoments<Volume>;
  using MomentsType = Volume;

  friend SeparatedMomentsCommon<SelfType,MomentsType>;

 public:
  using SeparatedMomentsCommon<SelfType, MomentsType>::SeparatedMomentsCommon;

  SeparatedMoments(void) = default;

  template <class GeometryType>
  static SeparatedMoments fillWithComplementMoments(
      const MomentsType& a_known_moments,
      const GeometryType& a_encompassing_geometry, const bool a_flipped);

  static SeparatedMoments fillWithComplementMoments(
      const MomentsType& a_known_moments,
      const MomentsType& a_encompassing_geometry_volume_moments,
      const bool a_flipped);

	/// \brief Set volumes from a list of 2 doubles.
	///
	/// Construct the volume moments from a list of 2 doubles.
	/// Expected order is Internal Volume, External Volume
	explicit SeparatedMoments(const double* a_list);

};

template <UnsignedIndex_t kArrayLength>
class SeparatedMoments<VolumeMomentsAndDoubles<kArrayLength>>
    : public SeparatedMomentsCommon<
          SeparatedMoments<VolumeMomentsAndDoubles<kArrayLength>>,
          VolumeMomentsAndDoubles<kArrayLength>> {
  using SelfType = SeparatedMoments<VolumeMomentsAndDoubles<kArrayLength>>;
  using MomentsType = VolumeMomentsAndDoubles<kArrayLength>;

  friend SeparatedMomentsCommon<SelfType,MomentsType>;

 public:
  using SeparatedMomentsCommon<SelfType, MomentsType>::SeparatedMomentsCommon;

  SeparatedMoments(void) = default;

  template <class GeometryType>
  static SeparatedMoments fillWithComplementMoments(
      const MomentsType& a_known_moments,
      const GeometryType& a_encompassing_geometry, const bool a_flipped);

  static SeparatedMoments fillWithComplementMoments(
      const MomentsType& a_known_moments,
      const MomentsType&
          a_encompassing_geometry_volume_moments,
      const bool a_flipped);

};

template <class MomentsType>
inline std::ostream& operator<<(
    std::ostream& out,
    const SeparatedMoments<MomentsType>& a_separated_volume_moments);

/// \brief Overload * operator to multiply the two geometric moments.
template <class MomentsType>
SeparatedMoments<MomentsType> operator*(
    const SeparatedMoments<MomentsType>& a_svm, const double a_multiplier);
/// \brief Overload * operator to multiply the two geometric moments.
template <class MomentsType>
SeparatedMoments<MomentsType> operator*(
    const double a_multiplier, const SeparatedMoments<MomentsType>& a_svm);

}  // namespace IRL

#include "src/moments/separated_volume_moments.tpp"

#endif  // SRC_MOMENTS_SEPARATED_VOLUME_MOMENTS_H_
