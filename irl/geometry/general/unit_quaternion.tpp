// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_UNIT_QUATERNION_TPP_
#define IRL_GEOMETRY_GENERAL_UNIT_QUATERNION_TPP_

namespace IRL {

template <class ScalarType>
inline UnitQuaternionBase<ScalarType>::UnitQuaternionBase(
    const ScalarType a_rotation_amount_in_radians,
    const NormalBase<ScalarType>& a_rotation_axis) {
  quat_m[0] = cos(a_rotation_amount_in_radians / static_cast<ScalarType>(2));
  ScalarType scaling =
      sin(a_rotation_amount_in_radians / static_cast<ScalarType>(2));
  quat_m[1] = scaling * a_rotation_axis[0];
  quat_m[2] = scaling * a_rotation_axis[1];
  quat_m[3] = scaling * a_rotation_axis[2];
}

template <class ScalarType>
inline UnitQuaternionBase<ScalarType>
UnitQuaternionBase<ScalarType>::fromFourElements(ScalarType a_q0,
                                                 ScalarType a_q1,
                                                 ScalarType a_q2,
                                                 ScalarType a_q3) {
  return UnitQuaternion(a_q0, a_q1, a_q2, a_q3);
}

template <class ScalarType>
inline UnitQuaternionBase<ScalarType>
UnitQuaternionBase<ScalarType>::fromFourElementsNormalized(ScalarType a_q0,
                                                           ScalarType a_q1,
                                                           ScalarType a_q2,
                                                           ScalarType a_q3) {
  auto quaternion_to_normalize =
      UnitQuaternionBase<ScalarType>(a_q0, a_q1, a_q2, a_q3);
  quaternion_to_normalize.normalize();
  return quaternion_to_normalize;
}

template <class ScalarType>
inline const ScalarType& UnitQuaternionBase<ScalarType>::operator[](
    const UnsignedIndex_t a_elem) const {
  assert(a_elem < 4);
  return quat_m[a_elem];
}

template <class ScalarType>
inline void UnitQuaternionBase<ScalarType>::normalize(void) {
  ScalarType magnitude = this->magnitude();
  for (auto& element : quat_m) {
    element /= magnitude;
  }
}

template <class ScalarType>
inline ScalarType UnitQuaternionBase<ScalarType>::magnitude(void) const {
  return sqrt(quat_m[0] * quat_m[0] + quat_m[1] * quat_m[1] +
              quat_m[2] * quat_m[2] + quat_m[3] * quat_m[3]);
}

template <class ScalarType>
inline UnitQuaternionBase<ScalarType> UnitQuaternionBase<ScalarType>::inverse(
    void) const {
  return UnitQuaternionBase<ScalarType>(quat_m[0], -quat_m[1], -quat_m[2],
                                        -quat_m[3]);
}

template <class ScalarType>
inline UnitQuaternionBase<ScalarType> UnitQuaternionBase<ScalarType>::operator*(
    const UnitQuaternionBase<ScalarType>& a_unit_quaternion) const {
  UnitQuaternionBase<ScalarType> unit_quat_to_return;
  unit_quat_to_return[0] =
      a_unit_quaternion[0] * quat_m[0] - a_unit_quaternion[1] * quat_m[1] -
      a_unit_quaternion[2] * quat_m[2] - a_unit_quaternion[3] * quat_m[3];

  unit_quat_to_return[1] =
      a_unit_quaternion[0] * quat_m[1] + a_unit_quaternion[1] * quat_m[0] -
      a_unit_quaternion[2] * quat_m[3] + a_unit_quaternion[3] * quat_m[2];

  unit_quat_to_return[2] =
      a_unit_quaternion[0] * quat_m[2] + a_unit_quaternion[1] * quat_m[3] +
      a_unit_quaternion[2] * quat_m[0] - a_unit_quaternion[3] * quat_m[1];

  unit_quat_to_return[3] =
      a_unit_quaternion[0] * quat_m[3] - a_unit_quaternion[1] * quat_m[2] +
      a_unit_quaternion[2] * quat_m[1] + a_unit_quaternion[3] * quat_m[0];

  unit_quat_to_return.normalize();
  return unit_quat_to_return;
}

template <class ScalarType>
inline NormalBase<ScalarType> UnitQuaternionBase<ScalarType>::operator*(
    const NormalBase<ScalarType>& a_normal) const {
  UnitQuaternionBase<ScalarType> rotated_quat =
      (*this) *
      UnitQuaternionBase<ScalarType>(static_cast<ScalarType>(0), a_normal[0],
                                     a_normal[1], a_normal[2]) *
      this->inverse();
  NormalBase<ScalarType> rotated_normal =
      NormalBase<ScalarType>(rotated_quat[1], rotated_quat[2], rotated_quat[3]);
  rotated_normal.normalize();
  return rotated_normal;
}

template <class ScalarType>
inline ReferenceFrameBase<ScalarType> UnitQuaternionBase<ScalarType>::operator*(
    const ReferenceFrameBase<ScalarType>& a_reference_frame) const {
  return ReferenceFrameBase<ScalarType>((*this) * a_reference_frame[0],
                                        (*this) * a_reference_frame[1],
                                        (*this) * a_reference_frame[2]);
}

template <class ScalarType>
inline ScalarType& UnitQuaternionBase<ScalarType>::operator[](
    const UnsignedIndex_t a_elem) {
  assert(a_elem < 4);
  return quat_m[a_elem];
}

template <class ScalarType>
inline UnitQuaternionBase<ScalarType>::UnitQuaternionBase(ScalarType a_q0,
                                                          ScalarType a_q1,
                                                          ScalarType a_q2,
                                                          ScalarType a_q3)
    : quat_m{a_q0, a_q1, a_q2, a_q3} {}

}  // namespace IRL

#endif  // IRL_GEOMETRY_GENERAL_UNIT_QUATERNION_TPP_
