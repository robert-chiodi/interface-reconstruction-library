// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2020 Robert Chiodi  <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_H_
#define IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_H_

#include <math.h>
#include <cassert>

#include "irl/data_structures/stack_vector.h"
#include "irl/geometry/general/polynomial.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/helpers/helper.h"
#include "irl/paraboloid_reconstruction/aligned_paraboloid.h"

namespace IRL {

// Paraboloid in the form
// z = a + b*x + c*y + d*x^2 + e*x*y + f*y^2
class Paraboloid : public Polynomial<6> {
 public:
  using Polynomial<6>::Polynomial;

  // Generates equivalent AlignedParboloid of the form
  // Ax''^2 + By''^2 = z'', and returns
  // translation and rotation required to move into the
  // R'' = {x'', y'', z''} coordinate system.
  AlignedParaboloid generateAlignedParaboloid(
      Pt* a_translation, UnitQuaternion* a_rotation) const {
    assert(a_translation != nullptr);
    assert(a_rotation != nullptr);

    const double theta =
        0.5 * atan(this->e() / (safelyTiny(this->d() - this->f())));
    const double cos_t = cos(theta);
    const double sin_t = sin(theta);
    AlignedParaboloid aligned_paraboloid;
    const double A = this->d() * cos_t * cos_t + this->f() * sin_t * sin_t +
                     this->e() * cos_t * sin_t;
    const double B = this->f() * cos_t * cos_t + this->d() * sin_t * sin_t -
                     this->e() * cos_t * sin_t;
    aligned_paraboloid.a() = A;
    aligned_paraboloid.b() = B;
    aligned_paraboloid.c() = 1.0;

    // Translation to coordinate system R' where aligned paraboloid valid
    // Translation is R' = {x' = x + u, y' = y + v, z' = z + w}

    const double denominator =
        safelyTiny(4.0 * this->d() * this->f() - this->e() * this->e());
    const double u =
        (2.0 * this->b() * this->f() - this->c() * this->e()) / denominator;
    const double v =
        -(this->b() * this->e() - 2.0 * this->d() * this->c()) / denominator;
    const double w = -(this->a() + (-this->b() * this->b() * this->f() +
                                    this->b() * this->c() * this->e() -
                                    this->c() * this->c() * this->d()) /
                                       denominator);
    (*a_translation) = Pt(u, v, w);

    // Rotation from R' to R'' about z' axis
    // Represent as a unit quaternion
    (*a_rotation) = UnitQuaternion(theta, Normal(0.0, 0.0, 1.0));

    return aligned_paraboloid;
  }

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

 private:
};

inline Normal getParaboloidSurfaceNormal(const AlignedParaboloid& a_paraboloid,
                                         const Pt& a_pt) {
  return {2.0 * a_paraboloid.a() * a_pt[0], 2.0 * a_paraboloid.b() * a_pt[1],
          1.0};
};

// Returns solution to quadratic equation solve.
// The smallest solution will always be first.
inline StackVector<double, 2> solveQuadratic(const double a, const double b,
                                             const double c) {
  double discriminant = b * b - 4.0 * a * c;
  // By preventing discriminant = 0, we avoid cases with intersections tangent
  // to the paraboloid
  if (discriminant > 0.0) {
    if (a != 0.0) {
      discriminant = std::sqrt(discriminant);
      const double q = -0.5 * (b + std::copysign(discriminant, b));
      if (b == 0.0 && c == 0.0) {
        return StackVector<double, 2>({0.0, 0.0});
      } else if (q == 0.0) {
        return StackVector<double, 2>({0.0});
      }
      const double sol1 = q / a;
      const double sol2 = c / q;
      return sol1 < sol2 ? StackVector<double, 2>({sol1, sol2})
                         : StackVector<double, 2>({sol2, sol1});

    } else {
      return StackVector<double, 2>({-c / b});
    }
  }
  return StackVector<double, 2>();
};

inline Pt projectPtAlongLineOntoParaboloid(
    const AlignedParaboloid& a_paraboloid, const Normal& a_line,
    const Pt& a_starting_pt) {
  // a_line should be normalized before passing in to make
  // these checks make sense
  const double a = (a_paraboloid.a() * a_line[0] * a_line[0] +
                    a_paraboloid.b() * a_line[1] * a_line[1]);
  const double b =
      (a_line[2] + 2.0 * a_paraboloid.a() * a_starting_pt[0] * a_line[0] +
       2.0 * a_paraboloid.b() * a_starting_pt[1] * a_line[1]);
  const double c = (a_starting_pt[2] +
                    a_paraboloid.a() * a_starting_pt[0] * a_starting_pt[0] +
                    a_paraboloid.b() * a_starting_pt[1] * a_starting_pt[1]);
  // check if starting point is on paraboloid (then solution = 0)
  if (std::fabs(c) < 5.0 * DBL_EPSILON) {
    return a_starting_pt;
  } else {
    const auto solutions = solveQuadratic(a, b, c);
    if (solutions.size() == 0) {
      std::cout << a_line << a_starting_pt << std::endl;
    }
    assert(solutions.size() > 0);
    if (solutions.size() == 1) {
      assert(solutions[0] >= 0.0);
      return a_starting_pt + a_line * solutions[0];
    } else {
      assert(std::max(solutions[0], solutions[1]) >= 0.0);
      const double distance_along_line =
          solutions[0] > 0.0 && solutions[1] > 0.0
              ? std::min(solutions[0], solutions[1])
              : std::max(solutions[0], solutions[1]);
      return a_starting_pt + a_line * distance_along_line;
    }
  }
};
}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_H_
