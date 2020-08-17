// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_VOLUME_MOMENTS_AND_NORMAL_H_
#define SRC_MOMENTS_VOLUME_MOMENTS_AND_NORMAL_H_

#include "src/geometry/general/normal.h"
#include "src/moments/volume_moments.h"
#include "src/parameters/defined_types.h"

namespace IRL {
/// \brief Class that contains VolumeMoments for a triangle and the triangle's
/// normal.
class VolumeMomentsAndNormal {
 public:
  /// \brief Default constructor.
  VolumeMomentsAndNormal(void) = default;

  /// \brief Constructor given VolumeMoments and a normal.
  VolumeMomentsAndNormal(const VolumeMoments& a_volume_moments,
                         const Normal& a_normal);

  /// \brief Obtain un-normalized VolumeMoments from the supplied geometry.
  template <class GeometryType>
  static VolumeMomentsAndNormal calculateMoments(GeometryType* a_geometry);

  static VolumeMomentsAndNormal fromScalarConstant(const double a_constant);

  /// \brief Return address of volume_moments_m.
  VolumeMoments& volumeMoments(void);

  /// \brief Return const address of volume_moments_m.
  const VolumeMoments& volumeMoments(void) const;

  /// \brief Return normal.
  Normal& normal(void);

  /// \brief Return const normal.
  const Normal& normal(void) const;

  /// \brief Normalize centroid by corresponding volume.
  void normalizeByVolume(void);

  /// \brief Multiply centroid by corresponding volume.
  void multiplyByVolume(void);

  /// \brief Overload += operator to update moments.
  VolumeMomentsAndNormal& operator+=(const VolumeMomentsAndNormal& a_rhs);

  /// \brief Overload *= operator to update moments.
  VolumeMomentsAndNormal& operator*=(const double& a_rhs);

  /// \brief Overload /= operator to update moments.
  VolumeMomentsAndNormal& operator/=(const double& a_rhs);

  /// \brief Overload assignment to assign constant value to moments.
  VolumeMomentsAndNormal& operator=(const double a_value);

  /// \brief Default destructor.
  ~VolumeMomentsAndNormal(void) = default;

 private:
  explicit VolumeMomentsAndNormal(const double a_constant);

  VolumeMoments volume_moments_m;
  Normal normal_m;
};

inline std::ostream& operator<<(
    std::ostream& out,
    const VolumeMomentsAndNormal& a_volume_moments_and_normal);

}  // namespace IRL

#include "src/moments/volume_moments_and_normal.tpp"

#endif  // SRC_MOMENTS_VOLUME_MOMENTS_AND_NORMAL_H_
