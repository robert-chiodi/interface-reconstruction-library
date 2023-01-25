// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_VOLUME_MOMENTS_H_
#define IRL_MOMENTS_VOLUME_MOMENTS_H_

#include <utility>

#include "irl/geometry/general/pt.h"
#include "irl/moments/volume.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

/// \brief Zeroeth (volume) and first order (centroid) geometric moments
template <class ScalarType>
class VolumeMomentsBase {
 public:
  using gradient_type = VolumeMomentsBase;
  using value_type = ScalarType;

  /// \brief Default constructor.
  VolumeMomentsBase(void);

  /// \brief Constructor that initializes volume and centroid.
  constexpr VolumeMomentsBase(const ScalarType a_volume,
                              const PtBase<ScalarType>& a_centroid);
  constexpr VolumeMomentsBase(const VolumeMomentsBase<double>& a_moments);
  constexpr VolumeMomentsBase(const VolumeMomentsBase<Quad_t>& a_moments);

  static constexpr VolumeMomentsBase fromRawDoublePointer(
      const ScalarType* a_list);

  static constexpr VolumeMomentsBase fromScalarConstant(
      const ScalarType a_value);

  /// \brief Obtain un-normalized VolumeMoments from the supplied geometry.
  template <class GeometryType>
  static VolumeMomentsBase calculateMoments(GeometryType* a_geometry);

  /// \brief Return value of stored volume.
  VolumeBase<ScalarType>& volume(void);
  /// \brief Return const reference to stored volume.
  const VolumeBase<ScalarType>& volume(void) const;
  /// \brief Return copy of stored centroid.
  PtBase<ScalarType>& centroid(void);
  /// \brief Return const reference to stored centroid.
  const PtBase<ScalarType>& centroid(void) const;

  /// \brief Divide the centroid by the volume.
  void normalizeByVolume(void);

  /// \brief Multiply the centroid by the volume.
  void multiplyByVolume(void);

  /// \brief Overload += operator to update moments.
  VolumeMomentsBase& operator+=(const VolumeMomentsBase& a_rhs);

  /// \brief Overload *= operator to multiply by constant double
  VolumeMomentsBase& operator*=(const ScalarType a_rhs);

  /// \brief Overload *= operator to multiply by constant double
  VolumeMomentsBase& operator/=(const ScalarType a_rhs);

  /// \brief Overload assignment to assign constant value to moments.
  VolumeMomentsBase& operator=(const ScalarType a_value);

  // Unary - operator
  VolumeMomentsBase operator-(void) const;

  /// \brief Default destructor.
  ~VolumeMomentsBase(void) = default;

 private:
  /// \brief Construct VolumeMoments from a list of doubles.
  ///
  /// Construct VolumeMoments from a list of 4 doubles. The necessary order
  /// is volume, centroid_x, centroid_y, centroid_z.
  constexpr explicit VolumeMomentsBase(const ScalarType* a_list);

  /// \brief Construct that initializes volume/centroid as a value.
  constexpr explicit VolumeMomentsBase(const ScalarType a_value);

  VolumeBase<ScalarType> volume_m;  ///< \brief Zeroeth moment (volume).
  PtBase<ScalarType> centroid_m;    ///< \brief First moment (centroid).
};

template <class ScalarType>
inline std::ostream& operator<<(
    std::ostream& out, const VolumeMomentsBase<ScalarType>& a_volume_moments);

/// \brief Overload + operator to add two geometric moments together
template <class ScalarType>
inline VolumeMomentsBase<ScalarType> operator+(
    const VolumeMomentsBase<ScalarType>& a_vm1,
    const VolumeMomentsBase<ScalarType>& a_vm2);
/// \brief Overload - operator to subtract one
/// geometric moment object from another.
template <class ScalarType>
inline VolumeMomentsBase<ScalarType> operator-(
    const VolumeMomentsBase<ScalarType>& a_vm1,
    const VolumeMomentsBase<ScalarType>& a_vm2);
/// \brief Overload * operator to multiply volume/centroid
template <class ScalarType>
inline VolumeMomentsBase<ScalarType> operator*(
    const ScalarType a_multiplier, const VolumeMomentsBase<ScalarType>& a_vm);
/// \brief Overload * operator to multiply volume/centroid
template <class ScalarType>
inline VolumeMomentsBase<ScalarType> operator*(
    const VolumeMomentsBase<ScalarType>& a_vm, const ScalarType a_multiplier);

using VolumeMoments = VolumeMomentsBase<double>;

}  // namespace IRL

#include "irl/moments/volume_moments.tpp"

#endif  // IRL_MOMENTS_VOLUME_MOMENTS_H_
