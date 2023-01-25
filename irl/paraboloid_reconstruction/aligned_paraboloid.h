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
#include <ostream>

#include "irl/geometry/general/plane.h"
#include "irl/geometry/general/polynomial.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/paraboloid_reconstruction/ellipse.h"

namespace IRL {

// Paraboloid in the form z + a*x^2 + b*y^2 = 0
template <class ScalarType>
class AlignedParaboloidBase : public PolynomialBase<2, ScalarType> {
 public:
  using PolynomialBase<2, ScalarType>::PolynomialBase;
  using value_type = ScalarType;

  constexpr AlignedParaboloidBase(
      const AlignedParaboloidBase<double>& a_aligned_paraboloid)
      : PolynomialBase<2, ScalarType>(std::array<ScalarType, 2>{
            static_cast<ScalarType>(a_aligned_paraboloid.a()),
            static_cast<ScalarType>(a_aligned_paraboloid.b())}){};
  constexpr AlignedParaboloidBase(
      const AlignedParaboloidBase<Quad_t>& a_aligned_paraboloid)
      : PolynomialBase<2, ScalarType>(std::array<ScalarType, 2>{
            static_cast<ScalarType>(a_aligned_paraboloid.a()),
            static_cast<ScalarType>(a_aligned_paraboloid.b())}){};

  AlignedParaboloidBase flipParaboloid(void) const {
    return AlignedParaboloidBase(
        std::array<ScalarType, 2>{{-this->a(), -this->b()}});
  }

  ScalarType& a(void) { return (*this)[0]; }
  ScalarType a(void) const { return (*this)[0]; }
  ScalarType& b(void) { return (*this)[1]; }
  ScalarType b(void) const { return (*this)[1]; }

  // Ellipse in the form Ax^2 + By^2 + Cxy + Dx + Ey + F = 0
  EllipseBase<ScalarType> intersectWithPlane(
      const PlaneBase<ScalarType>& a_plane) const {
    EllipseBase<ScalarType> ellipse_to_return;
    if constexpr (std::is_same<ScalarType, Quad_t>::value) {
      if (fabsq(a_plane.normal()[2]) < FLT128_EPSILON) {
        // no intersection
        return ellipse_to_return;
      }

    } else {
      if (std::fabs(a_plane.normal()[2]) < DBL_EPSILON) {
        // no intersection
        return ellipse_to_return;
      }
    }
    // IRL Plane defined as n_x x + n_y y + n_z z - d = 0
    ellipse_to_return.a() = this->a();
    ellipse_to_return.b() = this->b();
    ellipse_to_return.c() =
        static_cast<ScalarType>(0);  // AlignedParaboloid is aligned along x-y
    ellipse_to_return.d() = -a_plane.normal()[0] / a_plane.normal()[2];
    ellipse_to_return.e() = -a_plane.normal()[1] / a_plane.normal()[2];
    ellipse_to_return.f() = a_plane.distance() / a_plane.normal()[2];
    return ellipse_to_return;
  }

 private:
};

using AlignedParaboloid = AlignedParaboloidBase<double>;

template <class ScalarType>
inline std::ostream& operator<<(
    std::ostream& out,
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid) {
  out << "0 = " << a_aligned_paraboloid.a() << " * x^2 + "
      << a_aligned_paraboloid.b() << " * y^2 + z";
  return out;
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTIONS_ALIGNED_PARABOLOID_H_
