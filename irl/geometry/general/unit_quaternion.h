// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_UNIT_QUATERNION_H_
#define IRL_GEOMETRY_GENERAL_UNIT_QUATERNION_H_

#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/reference_frame.h"

namespace IRL {

/// \brief Unit quaternion to be used to perform rotations.
template <class ScalarType>
class UnitQuaternionBase {
 public:
  using value_type = ScalarType;
  UnitQuaternionBase(void) = default;

  /// \brief Construct quaternion given rotation amount (in radians) and
  /// rotation axis
  UnitQuaternionBase(const ScalarType a_rotation_amount_in_radians,
                     const NormalBase<ScalarType>& a_rotation_axis);

  static UnitQuaternionBase fromFourElements(ScalarType a_q0, ScalarType a_q1,
                                             ScalarType a_q2, ScalarType a_q3);

  static UnitQuaternionBase fromFourElementsNormalized(ScalarType a_q0,
                                                       ScalarType a_q1,
                                                       ScalarType a_q2,
                                                       ScalarType a_q3);

  /// \brief Const access through `overload[]`
  const ScalarType& operator[](const UnsignedIndex_t a_elem) const;

  /// \brief Normalize the unit quaternion
  inline void normalize(void);

  /// \brief Return magnitude of the quaternion.
  ScalarType magnitude(void) const;

  /// \brief Return a copy of the inverse UnitQuaternionBase
  inline UnitQuaternionBase inverse(void) const;

  /// \brief Compile unit quaternions to perform successive rotations.
  inline UnitQuaternionBase operator*(
      const UnitQuaternionBase& a_unit_quaternion) const;

  /// \brief Rotate a normal
  inline NormalBase<ScalarType> operator*(
      const NormalBase<ScalarType>& a_normal) const;

  /// \brief Rotate a reference frame
  inline ReferenceFrameBase<ScalarType> operator*(
      const ReferenceFrameBase<ScalarType>& a_reference_frame) const;

  /// \brief Default destructor
  ~UnitQuaternionBase(void) = default;

 private:
  /// \brief Access through `overload[]`
  ScalarType& operator[](const UnsignedIndex_t a_elem);

  /// \brief Constructor given 4 doubles
  UnitQuaternionBase(ScalarType a_q0, ScalarType a_q1, ScalarType a_q2,
                     ScalarType a_q3);

  /// \brief Storage for the 4 elements of a quaternion
  std::array<ScalarType, 4> quat_m;
};

using UnitQuaternion = UnitQuaternionBase<double>;

}  // namespace IRL

#include "irl/geometry/general/unit_quaternion.tpp"

#endif  // IRL_GEOMETRY_GENERAL_UNIT_QUATERNION_H_
