// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi  <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION__ELLIPSE_H_
#define IRL_PARABOLOID_RECONSTRUCTION__ELLIPSE_H_

#include <math.h>
#include <cassert>

#include "irl/geometry/general/polynomial.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/unit_quaternion.h"

namespace IRL {

// FIXME: Rewrite equation and use to be more like Paraboloid in
// that coefficients start from lowest order and then go higher

// Ellipse in the form Ax^2 + By^2 + Cxy + Dx + Ey + F = 0
class Ellipse : public Polynomial<6> {
 public:
  using Polynomial<6>::Polynomial;
  double& a(void) { return (*this)[0]; }
  double a(void) const { return (*this)[0]; }
  double& b(void) { return (*this)[1]; }
  double b(void) const { return (*this)[1]; }
  double& c(void) { return (*this)[2]; }
  double c(void) const { return (*this)[2]; }
  double& d(void) { return (*this)[3]; }
  double d(void) const { return (*this)[3]; }
  double& e(void) { return (*this)[4]; }
  double e(void) const { return (*this)[4]; }
  double& f(void) { return (*this)[5]; }
  double f(void) const { return (*this)[5]; }

  bool isEllipseReal(void) const {
    return -this->f() + (0.25 / this->a()) * this->d() * this->d() +
               (0.25 / this->b()) * this->e() * this->e() <
           0.0;
  }

  std::array<double, 2> calculateCentroid(void) const {
    std::array<double, 2> centroid_2D;
    assert(std::fabs(this->c()) < DBL_EPSILON);
    centroid_2D[0] = -0.5 * this->d() / this->a();
    centroid_2D[1] = -0.5 * this->e() / this->b();
    return centroid_2D;
  }

 private:
};

inline std::ostream& operator<<(std::ostream& out, const Ellipse& a_ellipse) {
  out << "Ellipse equation: " << a_ellipse.a() << " *x^2 + " << a_ellipse.b()
      << " *y^2 + " << a_ellipse.c() << " *xy + " << a_ellipse.d() << " *x + "
      << a_ellipse.e() << " *y + " << a_ellipse.f() << std::endl;
  return out;
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTION__ELLIPSE_H_
