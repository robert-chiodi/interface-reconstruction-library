// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/general/rotations.h"

namespace IRL {

UnitQuaternion rotateNormalOntoNormal(const Normal& a_normal_original,
                                      const Normal& a_normal_destination) {
  double rotation_amount;
  Normal rotation_axis;
  return rotateNormalOntoNormal(a_normal_original, a_normal_destination,
                                &rotation_amount, &rotation_axis);
}

UnitQuaternion rotateNormalOntoNormal(const Normal& a_normal_original,
                                      const Normal& a_normal_destination,
                                      double* a_rotation_amount,
                                      Normal* a_rotation_axis) {
  double ndot = a_normal_original * a_normal_destination;

  // If same normal, return zero rotation unit quaternion
  if (ndot > 1.0 - 1.0e-13) {
    *a_rotation_amount = 0.0;
    *a_rotation_axis = Normal(0.0, 0.0, 1.0);
    return UnitQuaternion(*a_rotation_amount, *a_rotation_axis);
  }

  // Check if vectors are anti-parallel
  if (ndot < -1.0 + 1.0e-13) {
    *a_rotation_axis = crossProduct(Normal(0.0, 0.0, 1.0), a_normal_original);
    // Check if a_normal original was the z-axis
    if (magnitude(*a_rotation_axis) < 1.0e-8) {
      *a_rotation_axis = crossProduct(Normal(0.0, 1.0, 0.0), a_normal_original);
    }
    a_rotation_axis->normalize();
    *a_rotation_amount = M_PI;
    return UnitQuaternion(*a_rotation_amount, *a_rotation_axis);
  }

  // Take cross product of the two vectors
  *a_rotation_axis =
      crossProductNormalized(a_normal_original, a_normal_destination);
  *a_rotation_amount = std::acos(ndot);
  return UnitQuaternion(*a_rotation_amount, *a_rotation_axis);
}

ReferenceFrame getOrthonormalSystem(const Normal& a_normal) {
  ReferenceFrame xyz_frame(Normal(1.0, 0.0, 0.0), Normal(0.0, 1.0, 0.0),
                           Normal(0.0, 0.0, 1.0));
  UnitQuaternion rotation = rotateNormalOntoNormal(xyz_frame[2], a_normal);
  return rotation * xyz_frame;
}

Normal getSharedNormal(const Normal& a_normal_0, const Normal& a_normal_1) {
  double half_rotation_angle;
  Normal rotation_axis;
  return getSharedNormal(a_normal_0, a_normal_1, &half_rotation_angle,
                         &rotation_axis);
}

Normal getSharedNormal(const Normal& a_normal_0, const Normal& a_normal_1,
                       double* a_half_rotation_angle, Normal* a_rotation_axis) {
  UnitQuaternion full_rotation = rotateNormalOntoNormal(
      a_normal_0, a_normal_1, a_half_rotation_angle, a_rotation_axis);
  *a_half_rotation_angle *= 0.5;
  UnitQuaternion half_rotation(*a_half_rotation_angle, *a_rotation_axis);
  return half_rotation * a_normal_0;
}

}  // namespace IRL
