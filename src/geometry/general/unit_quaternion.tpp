// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_UNIT_QUATERNION_TPP_
#define SRC_GEOMETRY_GENERAL_UNIT_QUATERNION_TPP_

namespace IRL {

inline UnitQuaternion::UnitQuaternion(const double a_rotation_amount_in_radians,
                                      const Normal& a_rotation_axis) {
  quat_m[0] = std::cos(0.5 * a_rotation_amount_in_radians);
  double scaling = std::sin(0.5 * a_rotation_amount_in_radians);
  quat_m[1] = scaling * a_rotation_axis[0];
  quat_m[2] = scaling * a_rotation_axis[1];
  quat_m[3] = scaling * a_rotation_axis[2];
}

inline UnitQuaternion UnitQuaternion::fromFourElements(double a_q0, double a_q1,
                                                       double a_q2,
                                                       double a_q3) {
  return UnitQuaternion(a_q0, a_q1, a_q2, a_q3);
}

inline UnitQuaternion UnitQuaternion::fromFourElementsNormalized(double a_q0,
                                                                 double a_q1,
                                                                 double a_q2,
                                                                 double a_q3) {
  auto quaternion_to_normalize = UnitQuaternion(a_q0, a_q1, a_q2, a_q3);
  quaternion_to_normalize.normalize();
  return quaternion_to_normalize;
}

inline const double& UnitQuaternion::operator[](
    const UnsignedIndex_t a_elem) const {
  assert(a_elem < 4);
  return quat_m[a_elem];
}

inline void UnitQuaternion::normalize(void) {
  double magnitude = this->magnitude();
  for (auto& element : quat_m) {
    element /= magnitude;
  }
}

inline double UnitQuaternion::magnitude(void) const {
  return std::sqrt(quat_m[0] * quat_m[0] + quat_m[1] * quat_m[1] +
                   quat_m[2] * quat_m[2] + quat_m[3] * quat_m[3]);
}

inline UnitQuaternion UnitQuaternion::inverse(void) const {
  return UnitQuaternion(quat_m[0], -quat_m[1], -quat_m[2], -quat_m[3]);
}

inline UnitQuaternion UnitQuaternion::operator*(
    const UnitQuaternion& a_unit_quaternion) const {
  UnitQuaternion unit_quat_to_return;
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

inline Normal UnitQuaternion::operator*(const Normal& a_normal) const {
  UnitQuaternion rotated_quat =
      (*this) * UnitQuaternion(0.0, a_normal[0], a_normal[1], a_normal[2]) *
      this->inverse();
  Normal rotated_normal =
      Normal(rotated_quat[1], rotated_quat[2], rotated_quat[3]);
  rotated_normal.normalize();
  return rotated_normal;
}

inline ReferenceFrame UnitQuaternion::operator*(
    const ReferenceFrame& a_reference_frame) const {
  return ReferenceFrame((*this) * a_reference_frame[0],
                        (*this) * a_reference_frame[1],
                        (*this) * a_reference_frame[2]);
}

inline double& UnitQuaternion::operator[](const UnsignedIndex_t a_elem) {
  assert(a_elem < 4);
  return quat_m[a_elem];
}

inline UnitQuaternion::UnitQuaternion(double a_q0, double a_q1, double a_q2,
                                      double a_q3)
    : quat_m{a_q0, a_q1, a_q2, a_q3} {}

}  // namespace IRL

#endif  // SRC_GEOMETRY_GENERAL_UNIT_QUATERNION_TPP_
