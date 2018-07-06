// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_UNIT_QUATERNION_H_
#define SRC_GEOMETRY_GENERAL_UNIT_QUATERNION_H_

#include "src/geometry/general/normal.h"
#include "src/geometry/general/reference_frame.h"

namespace IRL {

/// \brief Unit quaternion to be used to perform rotations.
class UnitQuaternion {
 public:
  UnitQuaternion(void) = default;

  /// \brief Construct quaternion given rotation amount (in radians) and
  /// rotation axis
  UnitQuaternion(const double a_rotation_amount_in_radians,
                 const Normal& a_rotation_axis);

  static UnitQuaternion fromFourElements(double a_q0, double a_q1, double a_q2,
                                         double a_q3);

  static UnitQuaternion fromFourElementsNormalized(double a_q0, double a_q1,
                                                   double a_q2, double a_q3);

  /// \brief Const access through `overload[]`
  const double& operator[](const UnsignedIndex_t a_elem) const;

  /// \brief Normalize the unit quaternion
  inline void normalize(void);

  /// \brief Return magnitude of the quaternion.
  double magnitude(void) const;

  /// \brief Return a copy of the inverse UnitQuaternion
  inline UnitQuaternion inverse(void) const;

  /// \brief Compile unit quaternions to perform successive rotations.
  inline UnitQuaternion operator*(
      const UnitQuaternion& a_unit_quaternion) const;

  /// \brief Rotate a normal
  inline Normal operator*(const Normal& a_normal) const;

  /// \brief Rotate a reference frame
  inline ReferenceFrame operator*(
      const ReferenceFrame& a_reference_frame) const;

  /// \brief Default destructor
  ~UnitQuaternion(void) = default;

 private:
  /// \brief Access through `overload[]`
  double& operator[](const UnsignedIndex_t a_elem);

  /// \brief Constructor given 4 doubles
  UnitQuaternion(double a_q0, double a_q1, double a_q2, double a_q3);

  /// \brief Storage for the 4 elements of a quaternion
  std::array<double, 4> quat_m;
};

}  // namespace IRL

#include "src/geometry/general/unit_quaternion.tpp"

#endif  // SRC_GEOMETRY_GENERAL_UNIT_QUATERNION_H_
