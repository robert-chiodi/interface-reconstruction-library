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
    const Normal& a_end_tangent, const Normal& a_plane_normal,
    const AlignedParaboloid& a_paraboloid) {
  // Set what can already be set
  start_point_m = a_start_pt;
  end_point_m = a_end_pt;
  start_point_id_m = reinterpret_cast<std::uintptr_t>(&a_start_pt);
  end_point_id_m = reinterpret_cast<std::uintptr_t>(&a_end_pt);

  // Compute control point
  const Normal edge_vector = end_point_m - start_point_m;
  const Pt average_pt = 0.5 * (start_point_m + end_point_m);
  const Normal n_cross_t0 = crossProduct(a_plane_normal, a_start_tangent);
  if (std::fabs(n_cross_t0 * a_end_tangent) < 100.0 * DBL_EPSILON) {
    control_point_m = average_pt;
    weight_m = DBL_MAX;
  } else {
    assert(std::fabs(n_cross_t0 * a_end_tangent) >= 100.0 * DBL_EPSILON);
    const double lambda_1 =
        -(n_cross_t0 * edge_vector) / (n_cross_t0 * a_end_tangent);
    control_point_m = Pt(end_point_m + lambda_1 * a_end_tangent);
    const double ct_correction =
        Normal(control_point_m - start_point_m) * a_plane_normal;
    control_point_m = control_point_m - ct_correction * a_plane_normal;
    auto mid_to_control = Normal(control_point_m - average_pt);
    if (squaredMagnitude(mid_to_control) < 1.0e4 * DBL_EPSILON * DBL_EPSILON) {
      weight_m = DBL_MAX;
    } else {
      mid_to_control.normalize();
      const double mid_dot_face_normal = mid_to_control * a_plane_normal;
      mid_to_control = mid_to_control - mid_dot_face_normal * a_plane_normal;
      const auto projected_pt = projectPtAlongHalfLineOntoParaboloid(
          a_paraboloid, mid_to_control, average_pt);
      weight_m = 1.0;
      if (squaredMagnitude(projected_pt - control_point_m) <
          1.0e4 * DBL_EPSILON * DBL_EPSILON) {
        weight_m = DBL_MAX;
      } else if (projected_pt[0] == DBL_MAX) {
        if (std::fabs(a_paraboloid.a() * average_pt[0] * average_pt[0] +
                      a_paraboloid.b() * average_pt[1] * average_pt[1] +
                      average_pt[2]) < 100.0 * DBL_EPSILON) {
          weight_m = 0.0;
        } else {
          weight_m = -1.0e15;
        }
      } else {
        double previous_best = 0.0;
        double sign = 1.0;
        for (UnsignedIndex_t d = 0; d < 3; ++d) {
          const double denominator = projected_pt[d] - control_point_m[d];
          if (std::fabs(denominator) > previous_best) {
            previous_best = std::fabs(denominator);
            weight_m =
                (start_point_m[d] + end_point_m[d] - 2.0 * projected_pt[d]) /
                (2.0 * denominator);
          }
        }
        if (previous_best < DBL_EPSILON) {
          weight_m = DBL_MAX;
        }
        weight_m = std::max(weight_m, 0.0);
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

template <class PtTypeWithGradient>
inline RationalBezierArcWithGradient<
    PtTypeWithGradient>::RationalBezierArcWithGradient(void) {
  start_point_m = PtTypeWithGradient();
  control_point_m = PtTypeWithGradient();
  end_point_m = PtTypeWithGradient();
  weight_m = 0.0;
  weight_gradient_m = PtTypeWithGradient::gradient_type();
}

template <class PtTypeWithGradient>
inline RationalBezierArcWithGradient<PtTypeWithGradient>::
    RationalBezierArcWithGradient(const PtTypeWithGradient& a_start_pt,
                                  const PtTypeWithGradient& a_start_tangent,
                                  const PtTypeWithGradient& a_end_pt,
                                  const PtTypeWithGradient& a_end_tangent,
                                  const PtTypeWithGradient& a_plane_normal,
                                  const AlignedParaboloid& a_paraboloid) {
  using gradient_type = typename PtTypeWithGradient::gradient_type;
  start_point_m = a_start_pt;
  end_point_m = a_end_pt;
  const Pt& start_point = a_start_pt.getPt();
  const auto& start_point_grad = a_start_pt.getData();
  const Pt& start_tangent = a_start_tangent.getPt();
  const auto& start_tangent_grad = a_start_tangent.getData();
  const Pt& end_point = a_end_pt.getPt();
  const auto& end_point_grad = a_end_pt.getData();
  const Pt& end_tangent = a_end_tangent.getPt();
  const auto& end_tangent_grad = a_end_tangent.getData();

  // Compute control point
  const auto edge_vector_withgrad = end_point_m - start_point_m;
  const Normal edge_vector = Normal::fromPt(edge_vector_withgrad.getPt());
  const auto& edge_vector_grad = edge_vector_withgrad.getData();
  const auto average_pt_withgrad = 0.5 * (start_point_m + end_point_m);
  const auto& average_pt = average_pt_withgrad.getPt();
  const auto& average_pt_grad = average_pt_withgrad.getData();
  const auto& plane_normal = a_plane_normal.getPt();
  const auto& plane_normal_grad = a_plane_normal.getData();
  const Pt n_cross_t0 = Pt(
      plane_normal[1] * start_tangent[2] - plane_normal[2] * start_tangent[1],
      plane_normal[2] * start_tangent[0] - plane_normal[0] * start_tangent[2],
      plane_normal[0] * start_tangent[1] - plane_normal[1] * start_tangent[0]);
  PtTypeWithGradient n_cross_t0_withgrad = PtTypeWithGradient(n_cross_t0);
  auto& n_cross_t0_grad = n_cross_t0_withgrad.getData();
  n_cross_t0_grad[0] = plane_normal_grad[1] * start_tangent[2] +
                       plane_normal[1] * start_tangent_grad[2] -
                       plane_normal_grad[2] * start_tangent[1] -
                       plane_normal[2] * start_tangent_grad[1];
  n_cross_t0_grad[1] = plane_normal_grad[2] * start_tangent[0] +
                       plane_normal[2] * start_tangent_grad[0] -
                       plane_normal_grad[0] * start_tangent[2] -
                       plane_normal[0] * start_tangent_grad[2];
  n_cross_t0_grad[2] = plane_normal_grad[0] * start_tangent[1] +
                       plane_normal[0] * start_tangent_grad[1] -
                       plane_normal_grad[1] * start_tangent[0] -
                       plane_normal[1] * start_tangent_grad[0];
  // assert(std::fabs(Normal::fromPt(n_cross_t0) * Normal::fromPt(end_tangent))
  // >
  //        DBL_EPSILON);
  const double n_cross_t0_dot_edge_vector = n_cross_t0[0] * edge_vector[0] +
                                            n_cross_t0[1] * edge_vector[1] +
                                            n_cross_t0[2] * edge_vector[2];
  const gradient_type n_cross_t0_dot_edge_vector_grad =
      n_cross_t0_grad[0] * edge_vector[0] +
      n_cross_t0[0] * edge_vector_grad[0] +
      n_cross_t0_grad[1] * edge_vector[1] +
      n_cross_t0[1] * edge_vector_grad[1] +
      n_cross_t0_grad[2] * edge_vector[2] + n_cross_t0[2] * edge_vector_grad[2];
  const double n_cross_t0_dot_end_tangent = n_cross_t0[0] * end_tangent[0] +
                                            n_cross_t0[1] * end_tangent[1] +
                                            n_cross_t0[2] * end_tangent[2];
  const gradient_type n_cross_t0_dot_end_tangent_grad =
      n_cross_t0_grad[0] * end_tangent[0] +
      n_cross_t0[0] * end_tangent_grad[0] +
      n_cross_t0_grad[1] * end_tangent[1] +
      n_cross_t0[1] * end_tangent_grad[1] +
      n_cross_t0_grad[2] * end_tangent[2] + n_cross_t0[2] * end_tangent_grad[2];
  const double lambda_1 =
      -(n_cross_t0_dot_edge_vector) / safelyEpsilon(n_cross_t0_dot_end_tangent);
  const gradient_type lambda_1_grad =
      -(n_cross_t0_dot_edge_vector_grad * n_cross_t0_dot_end_tangent -
        n_cross_t0_dot_edge_vector * n_cross_t0_dot_end_tangent_grad) /
      safelyEpsilon(n_cross_t0_dot_end_tangent * n_cross_t0_dot_end_tangent);
  const Pt control_point = end_point + lambda_1 * end_tangent;
  control_point_m = PtTypeWithGradient(control_point);
  auto& control_point_grad = control_point_m.getData();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    control_point_grad[d] = end_point_grad[d] + lambda_1 * end_tangent_grad[d] +
                            lambda_1_grad * end_tangent[d];
  }
  const auto mid_to_control_withgrad = control_point_m - average_pt_withgrad;
  const auto projected_pt_withgrad =
      projectPtAlongHalfLineOntoParaboloidWithGradient<PtTypeWithGradient>(
          a_paraboloid, mid_to_control_withgrad, average_pt_withgrad);
  const auto& projected_pt = projected_pt_withgrad.getPt();
  const auto& projected_pt_grad = projected_pt_withgrad.getData();
  weight_m = 1.0;
  if (squaredMagnitude(projected_pt - control_point) <
      DBL_EPSILON * DBL_EPSILON) {
    weight_m = DBL_MAX;
    weight_gradient_m = gradient_type();
  } else {
    double previous_best = -DBL_MAX;
    UnsignedIndex_t chosen_id = 0;
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      const double denominator = projected_pt[d] - control_point[d];
      if (std::fabs(denominator) > previous_best) {
        chosen_id = d;
        previous_best = std::fabs(denominator);
        weight_m = (start_point[d] + end_point[d] - 2.0 * projected_pt[d]) /
                   safelyEpsilon(2.0 * denominator);
      }
    }
    const double denominator =
        projected_pt[chosen_id] - control_point[chosen_id];
    const gradient_type den_grad =
        projected_pt_grad[chosen_id] - control_point_grad[chosen_id];
    weight_gradient_m =
        ((start_point_grad[chosen_id] + end_point_grad[chosen_id] -
          2.0 * projected_pt_grad[chosen_id]) *
             denominator -
         (start_point[chosen_id] + end_point[chosen_id] -
          2.0 * projected_pt[chosen_id]) *
             den_grad) /
        safelyEpsilon(2.0 * denominator * denominator);
  }
}

template <class PtTypeWithGradient>
inline RationalBezierArc RationalBezierArcWithGradient<PtTypeWithGradient>::arc(
    void) const {
  return RationalBezierArc(start_point_m.getPt(), control_point_m.getPt(),
                           end_point_m.getPt(), weight_m);
}

template <class PtTypeWithGradient>
inline double RationalBezierArcWithGradient<PtTypeWithGradient>::weight(
    void) const {
  return weight_m;
}

template <class PtTypeWithGradient>
inline const typename PtTypeWithGradient::gradient_type&
RationalBezierArcWithGradient<PtTypeWithGradient>::weight_gradient(void) const {
  return weight_gradient_m;
}

template <class PtTypeWithGradient>
inline const PtTypeWithGradient&
RationalBezierArcWithGradient<PtTypeWithGradient>::start_point(void) const {
  return start_point_m;
}

template <class PtTypeWithGradient>
inline const PtTypeWithGradient&
RationalBezierArcWithGradient<PtTypeWithGradient>::control_point(void) const {
  return control_point_m;
}

template <class PtTypeWithGradient>
inline const PtTypeWithGradient&
RationalBezierArcWithGradient<PtTypeWithGradient>::end_point(void) const {
  return end_point_m;
}

template <class PtTypeWithGradient>
inline double& RationalBezierArcWithGradient<PtTypeWithGradient>::weight(void) {
  return weight_m;
}

template <class PtTypeWithGradient>
inline typename PtTypeWithGradient::gradient_type&
RationalBezierArcWithGradient<PtTypeWithGradient>::weight_gradient(void) {
  return weight_gradient_m;
}

template <class PtTypeWithGradient>
inline PtTypeWithGradient&
RationalBezierArcWithGradient<PtTypeWithGradient>::start_point(void) {
  return start_point_m;
}

template <class PtTypeWithGradient>
inline PtTypeWithGradient&
RationalBezierArcWithGradient<PtTypeWithGradient>::control_point(void) {
  return control_point_m;
}

template <class PtTypeWithGradient>
inline PtTypeWithGradient&
RationalBezierArcWithGradient<PtTypeWithGradient>::end_point(void) {
  return end_point_m;
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
