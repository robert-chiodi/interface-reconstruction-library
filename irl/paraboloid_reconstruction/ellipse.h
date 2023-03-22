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
template <class ScalarType>
class EllipseBase : public PolynomialBase<6, ScalarType> {
 public:
  using PolynomialBase<6, ScalarType>::PolynomialBase;
  using value_type = ScalarType;

  ScalarType& a(void) { return (*this)[0]; }
  ScalarType a(void) const { return (*this)[0]; }
  ScalarType& b(void) { return (*this)[1]; }
  ScalarType b(void) const { return (*this)[1]; }
  ScalarType& c(void) { return (*this)[2]; }
  ScalarType c(void) const { return (*this)[2]; }
  ScalarType& d(void) { return (*this)[3]; }
  ScalarType d(void) const { return (*this)[3]; }
  ScalarType& e(void) { return (*this)[4]; }
  ScalarType e(void) const { return (*this)[4]; }
  ScalarType& f(void) { return (*this)[5]; }
  ScalarType f(void) const { return (*this)[5]; }

  bool isEllipseReal(void) const {
    return -this->f() +
               (static_cast<ScalarType>(1) /
                (static_cast<ScalarType>(4) * this->a())) *
                   this->d() * this->d() +
               (static_cast<ScalarType>(1) /
                (static_cast<ScalarType>(4) * this->b())) *
                   this->e() * this->e() <
           static_cast<ScalarType>(0);
  }

  std::array<ScalarType, 2> calculateCentroid(void) const {
    std::array<ScalarType, 2> centroid_2D;
    assert(fabs(this->c()) < DBL_EPSILON);
    centroid_2D[0] = -this->d() / (static_cast<ScalarType>(2) * this->a());
    centroid_2D[1] = -this->e() / (static_cast<ScalarType>(2) * this->b());
    return centroid_2D;
  }

 private:
};

using Ellipse = EllipseBase<double>;

template <class ScalarType>
inline std::ostream& operator<<(std::ostream& out,
                                const EllipseBase<ScalarType>& a_ellipse) {
  out << "Ellipse equation: " << a_ellipse.a() << " *x^2 + " << a_ellipse.b()
      << " *y^2 + " << a_ellipse.c() << " *xy + " << a_ellipse.d() << " *x + "
      << a_ellipse.e() << " *y + " << a_ellipse.f() << std::endl;
  return out;
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTION__ELLIPSE_H_
