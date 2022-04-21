// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_RATIONAL_BEZIER_ARC_TPP_
#define IRL_PARABOLOID_RECONSTRUCTION_RATIONAL_BEZIER_ARC_TPP_

namespace IRL {

inline RationalBezierArc::RationalBezierArc(void)
    : weight_m{0.0},
      start_point_m{0.0, 0.0, 0.0},
      start_point_id_m{0},
      control_point_m{0.0, 0.0, 0.0},
      end_point_m{0.0, 0.0, 0.0},
      end_point_id_m{0} {}

inline RationalBezierArc::RationalBezierArc(const Pt& a_start_pt,
                                            const Pt& a_control_pt,
                                            const Pt& a_end_pt,
                                            const double a_weight)
    : weight_m(a_weight),
      start_point_m(a_start_pt),
      start_point_id_m(reinterpret_cast<std::uintptr_t>(&a_start_pt)),
      control_point_m(a_control_pt),
      end_point_m(a_end_pt),
      end_point_id_m(reinterpret_cast<std::uintptr_t>(&a_end_pt)) {}

inline RationalBezierArc::RationalBezierArc(
    const Pt& a_start_pt, const Normal& a_start_tangent, const Pt& a_end_pt,
    const Normal& a_end_tangent, const Plane& a_plane,
    const AlignedParaboloid& a_paraboloid) {
  // Set what can already be set
  start_point_m = a_start_pt;
  end_point_m = a_end_pt;
  start_point_id_m = reinterpret_cast<std::uintptr_t>(&a_start_pt);
  end_point_id_m = reinterpret_cast<std::uintptr_t>(&a_end_pt);

  // Compute control point
  const Normal edge_vector = end_point_m - start_point_m;
  const Pt average_pt = 0.5 * (start_point_m + end_point_m);
  const Normal n_cross_t0 = crossProduct(a_plane.normal(), a_start_tangent);
  assert(std::fabs(n_cross_t0 * a_end_tangent) > DBL_EPSILON);
  const double lambda_1 =
      -(n_cross_t0 * edge_vector) / (n_cross_t0 * a_end_tangent);
  control_point_m = Pt(end_point_m + lambda_1 * a_end_tangent);
  const auto mid_to_control = Normal(control_point_m - average_pt);
  const auto projected_pt = projectPtAlongHalfLineOntoParaboloid(
      a_paraboloid, mid_to_control, average_pt);
  weight_m = 1.0;
  if (squaredMagnitude(projected_pt - control_point_m) <
      DBL_EPSILON * DBL_EPSILON) {
    weight_m = DBL_MAX;
  } else {
    double previous_best = -DBL_MAX;
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      const double denominator = projected_pt[d] - control_point_m[d];
      if (std::fabs(denominator) > previous_best) {
        previous_best = std::fabs(denominator);
        weight_m = (start_point_m[d] + end_point_m[d] - 2.0 * projected_pt[d]) /
                   (2.0 * denominator);
      }
    }
  }
}

inline const double& RationalBezierArc::weight(void) const { return weight_m; }

inline const Pt& RationalBezierArc::start_point(void) const {
  return start_point_m;
}

inline const Pt& RationalBezierArc::control_point(void) const {
  return control_point_m;
}

inline const Pt& RationalBezierArc::end_point(void) const {
  return end_point_m;
}

inline Pt RationalBezierArc::point(const double t) const {
  assert(t >= 0.0 && t <= 1.0);
  if (weight_m > 1.0e15) {
    if (t < 0.5) {
      return Pt((1.0 - 2.0 * t) * start_point_m + 2.0 * t * control_point_m);
    } else {
      return Pt((2.0 - 2.0 * t) * control_point_m +
                (2.0 * t - 1.0) * end_point_m);
    }
  } else {
    const double denominator =
        (1.0 - t) * (1.0 - t) + 2.0 * weight_m * t * (1.0 - t) + t * t;
    return Pt(((1.0 - t) * (1.0 - t) * start_point_m +
               2.0 * weight_m * t * (1.0 - t) * control_point_m +
               t * t * end_point_m) /
              denominator);
  }
}

inline Pt RationalBezierArc::derivative(const double t) const {
  assert(t >= 0.0 && t <= 1.0);
  if (weight_m > 1.0e15) {
    if (t < 0.5) {
      return Pt(2.0 * (control_point_m - start_point_m));
    } else {
      return Pt(2.0 * (end_point_m - control_point_m));
    }
  } else {
    double denominator =
        (1.0 - t) * (1.0 - t) + 2.0 * weight_m * t * (1.0 - t) + t * t;
    denominator *= denominator;
    return Pt((2.0 * (start_point_m - end_point_m) * (1.0 - weight_m) * t * t +
               4.0 * (start_point_m - control_point_m) * weight_m * t -
               2.0 * (start_point_m - end_point_m) * t -
               2.0 * weight_m * (start_point_m - control_point_m)) /
              denominator);
  }
}

inline double RationalBezierArc::arc_length(void) const {
  // 3-point quadrature rule for arc_length calculation
  const double t0 = 0.5 * (1.0 - sqrt(3.0 / 5.0));
  const double t1 = 0.5;
  const double t2 = 0.5 * (1.0 + sqrt(3.0 / 5.0));
  const double w0 = 5.0 / 18.0;
  const double w1 = 8.0 / 18.0;
  const double w2 = 5.0 / 18.0;
  Pt pt_0 = derivative(t0);
  Pt pt_1 = derivative(t1);
  Pt pt_2 = derivative(t2);
  const double norm0 =
      sqrt(pt_0[0] * pt_0[0] + pt_0[1] * pt_0[1] + pt_0[2] * pt_0[2]);
  const double norm1 =
      sqrt(pt_1[0] * pt_1[0] + pt_1[1] * pt_1[1] + pt_1[2] * pt_1[2]);
  const double norm2 =
      sqrt(pt_2[0] * pt_2[0] + pt_2[1] * pt_2[1] + pt_2[2] * pt_2[2]);
  return w0 * norm0 + w1 * norm1 + w2 * norm2;
}

inline std::uintptr_t RationalBezierArc::start_point_id(void) const {
  return start_point_id_m;
}

inline std::uintptr_t RationalBezierArc::end_point_id(void) const {
  return end_point_id_m;
}

inline RationalBezierArc& RationalBezierArc::operator+=(const Pt& a_rhs) {
  this->start_point() += a_rhs;
  this->control_point() += a_rhs;
  this->end_point() += a_rhs;
  return (*this);
}

inline RationalBezierArc RationalBezierArc::operator-(void) const {
  return RationalBezierArc(end_point_m, control_point_m, start_point_m,
                           weight_m);
}

inline double& RationalBezierArc::weight(void) { return weight_m; }

inline Pt& RationalBezierArc::start_point(void) { return start_point_m; }

inline Pt& RationalBezierArc::control_point(void) { return control_point_m; }

inline Pt& RationalBezierArc::end_point(void) { return end_point_m; }

inline std::uintptr_t RationalBezierArc::start_point_id(void) {
  return start_point_id_m;
}

inline std::uintptr_t RationalBezierArc::end_point_id(void) {
  return end_point_id_m;
}

inline std::ostream& operator<<(
    std::ostream& out, const RationalBezierArc& a_rational_bezier_arc) {
  out.precision(16);
  out << std::scientific << a_rational_bezier_arc.weight() << " ";
  out << std::scientific << a_rational_bezier_arc.start_point() << " ";
  out << std::scientific << a_rational_bezier_arc.control_point() << " ";
  out << std::scientific << a_rational_bezier_arc.end_point() << " ";
  out << a_rational_bezier_arc.start_point_id() << " ";
  out << a_rational_bezier_arc.end_point_id();
  return out;
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTION_RATIONAL_BEZIER_ARC_TPP_
