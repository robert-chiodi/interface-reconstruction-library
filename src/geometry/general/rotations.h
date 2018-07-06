// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_ROTATIONS_H_
#define SRC_GEOMETRY_GENERAL_ROTATIONS_H_

#include <cassert>
#include <cmath>

#include "src/geometry/general/normal.h"
#include "src/geometry/general/reference_frame.h"
#include "src/geometry/general/unit_quaternion.h"

namespace IRL {

/// \file rotations.h
/// This file contains classes and functions
/// that deal with getting reference frames
/// and performing projections and rotations
/// into new spaces.
///
/// First, class definitions will appear.
/// Below that, inlined function definitions will be given.

/// \brief Return the quaternion that would rotate `a_normal_original` onto
/// `a_normal_destination`.
///
/// This function returns a quaternion that will rotate the normal
/// `a_normal_original` onto the normal `a_normal_destination`. Additionally,
/// the rotation amount (in radians) and the rotation axis will be given through
/// the pointers `a_rotation_amount` and `a_rotation_axis`, respectively.
///
/// \param[in] a_normal_original Original normal we want to rotate.
/// \param[in] a_normal_destination Normal we want original normal to become
/// from rotation. \param[out] a_rotation_amount Amount of rotation used to make
/// the returned `UnitQuaternion`. \param[out] a_rotation_axis Rotation axis
/// used to make the returned `UnitQuaternion`.
UnitQuaternion rotateNormalOntoNormal(const Normal& a_normal_original,
                                      const Normal& a_normal_destination,
                                      double* a_rotation_amount,
                                      Normal* a_rotation_axis);

/// \brief Same as other function, however does not return the angle of rotation
/// or rotation axis.
UnitQuaternion rotateNormalOntoNormal(const Normal& a_normal_original,
                                      const Normal& a_normal_destination);

/// \brief Create an orthonormal reference frame with axis 3 being a_normal.
ReferenceFrame getOrthonormalSystem(const Normal& a_normal);

/// \brief Get normal that is halfway between `a_normal_0` and `a_normal_1`.
/// Return through pointer the rotation to the shared normal and the axis of
/// rotation
Normal getSharedNormal(const Normal& a_normal_0, const Normal& a_normal_1,
                       double* a_half_rotation_angle, Normal* a_rotation_axis);

/// \brief Same as other function, however does not return the half-angle of
/// rotation or rotation axis.
Normal getSharedNormal(const Normal& a_normal_0, const Normal& a_normal_1);

}  // namespace IRL

#endif  // SRC_GEOMETRY_GENERAL_ROTATIONS_H_
