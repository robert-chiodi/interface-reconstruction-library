// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_VOLUME_MOMENTS_H_
#define SRC_MOMENTS_VOLUME_MOMENTS_H_

#include <utility>

#include "src/geometry/general/pt.h"
#include "src/moments/volume.h"
#include "src/parameters/defined_types.h"

namespace IRL {

/// \brief Zeroeth (volume) and first order (centroid) geometric moments
class VolumeMoments {
 public:
  /// \brief Default constructor.
  VolumeMoments(void);

  /// \brief Constructor that initializes volume and centroid.
  constexpr VolumeMoments(const double a_volume, const Pt& a_centroid);

  static constexpr VolumeMoments fromRawDoublePointer(const double* a_list);

  static constexpr VolumeMoments fromScalarConstant(const double a_value);

  /// \brief Obtain un-normalized VolumeMoments from the supplied geometry.
  template <class GeometryType>
  static VolumeMoments calculateMoments(GeometryType* a_geometry);

  /// \brief Return value of stored volume.
  Volume& volume(void);
  /// \brief Return const reference to stored volume.
  const Volume& volume(void) const;
  /// \brief Return copy of stored centroid.
  Pt& centroid(void);
  /// \brief Return const reference to stored centroid.
  const Pt& centroid(void) const;

  /// \brief Divide the centroid by the volume.
  void normalizeByVolume(void);

  /// \brief Multiply the centroid by the volume.
  void multiplyByVolume(void);

  /// \brief Overload += operator to update moments.
  VolumeMoments& operator+=(const VolumeMoments& a_rhs);

  /// \brief Overload *= operator to multiply by constant double
  VolumeMoments& operator*=(const double a_rhs);

  /// \brief Overload *= operator to multiply by constant double
  VolumeMoments& operator/=(const double a_rhs);

  /// \brief Overload assignment to assign constant value to moments.
  VolumeMoments& operator=(const double a_value);

  /// \brief Default destructor.
  ~VolumeMoments(void) = default;

 private:
  /// \brief Construct VolumeMoments from a list of doubles.
  ///
  /// Construct VolumeMoments from a list of 4 doubles. The necessary order
  /// is volume, centroid_x, centroid_y, centroid_z.
  constexpr explicit VolumeMoments(const double* a_list);

  /// \brief Construct that initializes volume/centroid as a value.
  constexpr explicit VolumeMoments(const double a_value);

  Volume volume_m;  ///< \brief Zeroeth moment (volume).
  Pt centroid_m;    ///< \brief First moment (centroid).
};

inline std::ostream& operator<<(std::ostream& out,
                                const VolumeMoments& a_volume_moments);

/// \brief Overload + operator to add two geometric moments together
inline VolumeMoments operator+(const VolumeMoments& a_vm1,
                               const VolumeMoments& a_vm2);
/// \brief Overload - operator to subtract one
/// geometric moment object from another.
inline VolumeMoments operator-(const VolumeMoments& a_vm1,
                               const VolumeMoments& a_vm2);
/// \brief Overload * operator to multiply volume/centroid
inline VolumeMoments operator*(const double a_multiplier,
                               const VolumeMoments& a_vm);
/// \brief Overload * operator to multiply volume/centroid
inline VolumeMoments operator*(const VolumeMoments& a_vm,
                               const double a_multiplier);

}  // namespace IRL

#include "src/moments/volume_moments.tpp"

#endif  // SRC_MOMENTS_VOLUME_MOMENTS_H_
