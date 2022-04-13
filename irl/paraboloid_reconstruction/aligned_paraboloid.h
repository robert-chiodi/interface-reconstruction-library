// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi  <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTIONS_ALIGNED_PARABOLOID_H_
#define IRL_PARABOLOID_RECONSTRUCTIONS_ALIGNED_PARABOLOID_H_

#include <math.h>
#include <cassert>

#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/polynomial.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/paraboloid_reconstruction/ellipse.h"

namespace IRL {

// Paraboloid in the form z + c + a*x^2 + b*y^2 = 0
class AlignedParaboloid : public Polynomial<3> {
 public:
  using Polynomial<3>::Polynomial;

  AlignedParaboloid flipParaboloid(void) const {
    return AlignedParaboloid(
        std::array<double, 3>{{-this->a(), -this->b(), -this->c()}});
  }

  double& a(void) { return (*this)[0]; }
  double a(void) const { return (*this)[0]; }
  double& b(void) { return (*this)[1]; }
  double b(void) const { return (*this)[1]; }
  double& c(void) { return (*this)[2]; }
  double c(void) const { return (*this)[2]; }

  // Ellipse in the form Ax^2 + By^2 + Cxy + Dx + Ey + F = 0
  Ellipse intersectWithPlane(const Plane& a_plane) const {
    Ellipse ellipse_to_return;
    if (std::fabs(a_plane.normal()[2]) < DBL_EPSILON) {
      // no intersection
      return ellipse_to_return;
    }
    // IRL Plane defined as n_x x + n_y y + n_z z - d = 0
    ellipse_to_return.a() = this->a();
    ellipse_to_return.b() = this->b();
    ellipse_to_return.c() = 0.0;  // AlignedParaboloid is aligned along x-y
    ellipse_to_return.d() = -a_plane.normal()[0] / a_plane.normal()[2];
    ellipse_to_return.e() = -a_plane.normal()[1] / a_plane.normal()[2];
    ellipse_to_return.f() = a_plane.distance() / a_plane.normal()[2];
    return ellipse_to_return;
  }

 private:
};

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTIONS_ALIGNED_PARABOLOID_H_
