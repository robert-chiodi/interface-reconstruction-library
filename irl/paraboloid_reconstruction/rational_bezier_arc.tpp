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

template <class ScalarType>
inline RationalBezierArcBase<ScalarType>::RationalBezierArcBase(void)
    : weight_m{static_cast<ScalarType>(0)},
      start_point_m{static_cast<ScalarType>(0), static_cast<ScalarType>(0),
                    static_cast<ScalarType>(0)},
      start_point_id_m{static_cast<ScalarType>(0)},
      control_point_m{static_cast<ScalarType>(0), static_cast<ScalarType>(0),
                      static_cast<ScalarType>(0)},
      end_point_m{static_cast<ScalarType>(0), static_cast<ScalarType>(0),
                  static_cast<ScalarType>(0)},
      end_point_id_m{0} {}

template <class ScalarType>
inline RationalBezierArcBase<ScalarType>::RationalBezierArcBase(
    const PtBase<ScalarType>& a_start_pt,
    const PtBase<ScalarType>& a_control_pt, const PtBase<ScalarType>& a_end_pt,
    const ScalarType a_weight)
    : weight_m(a_weight),
      start_point_m(a_start_pt),
      start_point_id_m(reinterpret_cast<std::uintptr_t>(&a_start_pt)),
      control_point_m(a_control_pt),
      end_point_m(a_end_pt),
      end_point_id_m(reinterpret_cast<std::uintptr_t>(&a_end_pt)) {}

template <class ScalarType>
inline RationalBezierArcBase<ScalarType>::RationalBezierArcBase(
    const PtBase<ScalarType>& a_start_pt,
    const PtBase<ScalarType>& a_control_pt, const PtBase<ScalarType>& a_end_pt,
    const std::uintptr_t a_start_id, const std::uintptr_t a_end_id,
    const ScalarType a_weight)
    : weight_m(a_weight),
      start_point_m(a_start_pt),
      start_point_id_m(a_start_id),
      control_point_m(a_control_pt),
      end_point_m(a_end_pt),
      end_point_id_m(a_end_id) {}

template <class ScalarType>
inline RationalBezierArcBase<ScalarType>::RationalBezierArcBase(
    const PtBase<ScalarType>& a_start_pt,
    const NormalBase<ScalarType>& a_control_pt,
    const PtBase<ScalarType>& a_end_pt,
    const NormalBase<ScalarType>& a_plane_normal,
    const AlignedParaboloidBase<ScalarType>& a_paraboloid) {
  /* Defining constants and types */
  const ScalarType EPSILON = machine_epsilon<ScalarType>();
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType HALF = static_cast<ScalarType>(0.5);
  const ScalarType ONEHUNDRED = static_cast<ScalarType>(100);
  const ScalarType DISTANCE_EPSILON = ONEHUNDRED * EPSILON;

  /* Function */
  // Set what can already be set
  start_point_m = a_start_pt;
  control_point_m = a_control_pt;
  end_point_m = a_end_pt;
  start_point_id_m = reinterpret_cast<std::uintptr_t>(&a_start_pt);
  end_point_id_m = reinterpret_cast<std::uintptr_t>(&a_end_pt);

  /* By caculating the end-point curvature
This calculate the curvature of the conic with (Hartmann1996)
and uses it to compute the weight (Farin1992)*/
  const ScalarType L = squaredMagnitude(a_control_pt - a_start_pt);
  if (L < DISTANCE_EPSILON * DISTANCE_EPSILON) {
    weight_m = ZERO;
  } else {
    const NormalBase<ScalarType> start_normal =
        getParaboloidSurfaceNormal(a_paraboloid, a_start_pt);
    const NormalBase<ScalarType> start_cross_prod =
        crossProduct(a_plane_normal, start_normal);
    const ScalarType cross_sq_0 = start_cross_prod[0] * start_cross_prod[0];
    const ScalarType cross_sq_1 = start_cross_prod[1] * start_cross_prod[1];
    const ScalarType cross_sq_2 = start_cross_prod[2] * start_cross_prod[2];
    const ScalarType D =
        fabs(a_paraboloid.a() * cross_sq_0 + a_paraboloid.b() * cross_sq_1);
    if (D < DISTANCE_EPSILON * DISTANCE_EPSILON) {
      weight_m = ZERO;
    } else {
      const ScalarType R = (cross_sq_0 + cross_sq_1 + cross_sq_2) / L;
      const ScalarType A =
          sqrt(R * R * R *
               squaredMagnitude(crossProduct(a_end_pt - a_start_pt,
                                             a_control_pt - a_end_pt)));
      if (A < DISTANCE_EPSILON * DISTANCE_EPSILON) {
        weight_m = ZERO;
      } else {
        weight_m = HALF * sqrt(A / D);
      }
    }
  }
}

template <class ScalarType>
inline RationalBezierArcBase<ScalarType>::RationalBezierArcBase(
    const PtBase<ScalarType>& a_start_pt,
    const NormalBase<ScalarType>& a_start_tangent,
    const PtBase<ScalarType>& a_end_pt,
    const NormalBase<ScalarType>& a_end_tangent,
    const NormalBase<ScalarType>& a_plane_normal,
    const AlignedParaboloidBase<ScalarType>& a_paraboloid) {
  /* Defining constants and types */
  const ScalarType EPSILON = machine_epsilon<ScalarType>();
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType HALF = ONE / TWO;
  const ScalarType ONEANDHALF = ONE + HALF;
  const ScalarType ONEHUNDRED = static_cast<ScalarType>(100);
  const ScalarType DISTANCE_EPSILON = ONEHUNDRED * EPSILON;

  /* Function */
  // Set what can already be set
  start_point_m = a_start_pt;
  end_point_m = a_end_pt;
  start_point_id_m = reinterpret_cast<std::uintptr_t>(&a_start_pt);
  end_point_id_m = reinterpret_cast<std::uintptr_t>(&a_end_pt);

  // if constexpr (std::is_same_v<ScalarType, double>) {
  //   /* Compute control point using 3 plane equations */
  //   const NormalBase<ScalarType> start_normal =
  //       getParaboloidSurfaceNormal(a_paraboloid, a_start_pt);
  //   const NormalBase<ScalarType> end_normal =
  //       getParaboloidSurfaceNormal(a_paraboloid, a_end_pt);
  //   Eigen::Matrix3f A;
  //   Eigen::Vector3f b;
  //   A << start_normal[0], start_normal[1], start_normal[2], end_normal[0],
  //       end_normal[1], end_normal[2], a_plane_normal[0], a_plane_normal[1],
  //       a_plane_normal[2];
  //   b << start_normal * a_start_pt, end_normal * a_end_pt,
  //       a_plane_normal * a_start_pt;
  //   Eigen::Vector3f x = A.colPivHouseholderQr().solve(b);
  //   if (!b.isApprox(A * x, 100.0 * DBL_EPSILON)) {
  //     control_point_m = HALF * (start_point_m + end_point_m);
  //     weight_m = -static_cast<ScalarType>(DBL_MAX);
  //   } else {
  //     control_point_m = PtBase<ScalarType>(x(0), x(1), x(2));

  // double determinant =
  //     start_normal[0] *
  //         (end_normal[1] * a_plane_normal[2] - a_plane_normal[1]) -
  //     start_normal[1] *
  //         (end_normal[0] * a_plane_normal[2] - a_plane_normal[0]) +
  //     (end_normal[0] * a_plane_normal[1] - end_normal[1] *
  //     a_plane_normal[0]);
  // double b0 = start_normal * a_start_pt;
  // double b1 = end_normal * a_end_pt;
  // double b2 = a_plane_normal * a_start_pt;

  // if (fabs(determinant) < DBL_EPSILON) {
  //   control_point_m = HALF * (start_point_m + end_point_m);
  //   weight_m = -static_cast<ScalarType>(DBL_MAX);
  // } else {
  //   double invdet = ONE / determinant;
  //   double x =
  //       invdet *
  //       (b0 * (end_normal[1] * a_plane_normal[2] - a_plane_normal[1]) -
  //        b1 * (start_normal[1] * a_plane_normal[2] - a_plane_normal[1]) +
  //        b2 * (start_normal[1] - end_normal[1]));
  //   double y =
  //       invdet *
  //       (-b0 * (end_normal[0] * a_plane_normal[2] - a_plane_normal[0]) +
  //        b1 * (start_normal[0] * a_plane_normal[2] - a_plane_normal[0]) -
  //        b2 * (start_normal[0] - end_normal[0]));
  //   double z = invdet * (b0 * (end_normal[0] * a_plane_normal[1] -
  //                              a_plane_normal[0] * end_normal[1]) -
  //                        b1 * (start_normal[0] * a_plane_normal[1] -
  //                              a_plane_normal[0] * start_normal[1]) +
  //                        b2 * (start_normal[0] * end_normal[1] -
  //                              end_normal[0] * start_normal[1]));

  //     /* By caculating the end-point curvature
  //  This calculate the curvature of the conic with (Hartmann1996)
  //  and uses it to compute the weight (Farin1992)*/
  //     const NormalBase<ScalarType> start_normal =
  //         getParaboloidSurfaceNormal(a_paraboloid, a_start_pt);
  //     const NormalBase<ScalarType> start_cross_prod =
  //         crossProduct(a_plane_normal, start_normal);
  //     const ScalarType cross_sq_0 = start_cross_prod[0] *
  //     start_cross_prod[0]; const ScalarType cross_sq_1 = start_cross_prod[1]
  //     * start_cross_prod[1]; const ScalarType cross_sq_2 =
  //     start_cross_prod[2] * start_cross_prod[2]; const ScalarType D =
  //     safelyTiny(
  //         fabs(a_paraboloid.a() * cross_sq_0 + a_paraboloid.b() *
  //         cross_sq_1));
  //     const ScalarType R =
  //         (cross_sq_0 + cross_sq_1 + cross_sq_2) /
  //         safelyTiny(squaredMagnitude(control_point_m - a_start_pt));
  //     const ScalarType A = squaredMagnitude(
  //         crossProduct(a_end_pt - a_start_pt, control_point_m - a_end_pt));
  //     weight_m = HALF * sqrt(sqrt(A * R * R * R) / D);
  //   }
  // } else {
  /* Compute control point using 2 tangents */
  const NormalBase<ScalarType> edge_vector = end_point_m - start_point_m;
  const PtBase<ScalarType> average_pt = HALF * (start_point_m + end_point_m);
  const NormalBase<ScalarType> n_cross_t0 =
      crossProduct(a_plane_normal, a_start_tangent);
  if (fabs(n_cross_t0 * a_end_tangent) < ONEHUNDRED * EPSILON) {
    control_point_m = average_pt;
    weight_m = static_cast<ScalarType>(DBL_MAX);
  } else {
    assert(fabs(n_cross_t0 * a_end_tangent) >= ONEHUNDRED * EPSILON);
    const ScalarType lambda_1 =
        -(n_cross_t0 * edge_vector) / (n_cross_t0 * a_end_tangent);
    control_point_m =
        PtBase<ScalarType>(end_point_m + lambda_1 * a_end_tangent);
    const ScalarType ct_correction =
        NormalBase<ScalarType>(control_point_m - start_point_m) *
        a_plane_normal;
    control_point_m = control_point_m - ct_correction * a_plane_normal;

    /* By caculating the end-point curvature
 This calculate the curvature of the conic with (Hartmann1996)
 and uses it to compute the weight (Farin1992)*/
    const ScalarType L = squaredMagnitude(control_point_m - a_start_pt);
    if (L < DISTANCE_EPSILON * DISTANCE_EPSILON) {
      weight_m = static_cast<ScalarType>(0);
    } else {
      const NormalBase<ScalarType> start_normal =
          getParaboloidSurfaceNormal(a_paraboloid, a_start_pt);
      const NormalBase<ScalarType> start_cross_prod =
          crossProduct(a_plane_normal, start_normal);
      const ScalarType cross_sq_0 = start_cross_prod[0] * start_cross_prod[0];
      const ScalarType cross_sq_1 = start_cross_prod[1] * start_cross_prod[1];
      const ScalarType cross_sq_2 = start_cross_prod[2] * start_cross_prod[2];
      const ScalarType D =
          fabs(a_paraboloid.a() * cross_sq_0 + a_paraboloid.b() * cross_sq_1);
      if (D < DISTANCE_EPSILON * DISTANCE_EPSILON) {
        weight_m = static_cast<ScalarType>(0);
      } else {
        const ScalarType R = (cross_sq_0 + cross_sq_1 + cross_sq_2) / L;
        const ScalarType A =
            sqrt(R * R * R *
                 squaredMagnitude(crossProduct(a_end_pt - a_start_pt,
                                               control_point_m - a_end_pt)));
        if (A < DISTANCE_EPSILON * DISTANCE_EPSILON) {
          weight_m = static_cast<ScalarType>(0);
        } else {
          weight_m = HALF * sqrt(A / D);
        }
      }
    }

    // /* By caculating position on arc at t=1/2 */
    // auto mid_to_control = NormalBase<ScalarType>(control_point_m -
    // average_pt); if (squaredMagnitude(mid_to_control) <
    //     ONEHUNDRED * ONEHUNDRED * EPSILON * EPSILON) {
    //   weight_m = static_cast<ScalarType>(DBL_MAX);
    // } else {
    //   const ScalarType crude_invmag =
    //       ONE / (fabs(mid_to_control[0]) + fabs(mid_to_control[1]) +
    //              fabs(mid_to_control[2]));
    //   mid_to_control *= crude_invmag;
    //   mid_to_control.approximatelyNormalize();
    //   const ScalarType normal_correction = mid_to_control * a_plane_normal;
    //   mid_to_control = mid_to_control - normal_correction * a_plane_normal;
    //   const auto projected_pt =
    //       projectPtAlongHalfLineOntoParaboloid<ScalarType>(
    //           a_paraboloid, mid_to_control, average_pt);
    //   weight_m = ONE;
    //   if (squaredMagnitude(projected_pt - control_point_m) <
    //       ONEHUNDRED * ONEHUNDRED * EPSILON * EPSILON) {
    //     weight_m = static_cast<ScalarType>(DBL_MAX);
    //   } else if (projected_pt[0] == static_cast<ScalarType>(DBL_MAX)) {
    //     if (fabs(a_paraboloid.a() * average_pt[0] * average_pt[0] +
    //              a_paraboloid.b() * average_pt[1] * average_pt[1] +
    //              average_pt[2]) < ONEHUNDRED * EPSILON) {
    //       weight_m = ZERO;
    //     } else {
    //       weight_m = -static_cast<ScalarType>(DBL_MAX);
    //     }
    //   } else {
    //     ScalarType previous_best = ZERO;
    //     ScalarType sign = ONE;
    //     for (UnsignedIndex_t d = 0; d < 3; ++d) {
    //       const ScalarType denominator = projected_pt[d] -
    //       control_point_m[d]; if (fabs(denominator) > previous_best) {
    //         previous_best = fabs(denominator);
    //         weight_m =
    //             (start_point_m[d] + end_point_m[d] - TWO * projected_pt[d]) /
    //             (TWO * denominator);
    //       }
    //     }
    //     if (previous_best < EPSILON) {
    //       weight_m = static_cast<ScalarType>(DBL_MAX);
    //     }
    //     weight_m = maximum(weight_m, ZERO);
    //   }
    // }
    // }
  }
}

template <class ScalarType>
inline const ScalarType& RationalBezierArcBase<ScalarType>::weight(void) const {
  return weight_m;
}

template <class ScalarType>
inline const PtBase<ScalarType>& RationalBezierArcBase<ScalarType>::start_point(
    void) const {
  return start_point_m;
}

template <class ScalarType>
inline const PtBase<ScalarType>&
RationalBezierArcBase<ScalarType>::control_point(void) const {
  return control_point_m;
}

template <class ScalarType>
inline const PtBase<ScalarType>& RationalBezierArcBase<ScalarType>::end_point(
    void) const {
  return end_point_m;
}

template <class ScalarType>
inline PtBase<ScalarType> RationalBezierArcBase<ScalarType>::point(
    const ScalarType t) const {
  assert(t >= 0.0 && t <= 1.0);
  if (weight_m > 1.0e15) {
    if (t < 0.5) {
      return PtBase<ScalarType>((1.0 - 2.0 * t) * start_point_m +
                                2.0 * t * control_point_m);
    } else {
      return PtBase<ScalarType>((2.0 - 2.0 * t) * control_point_m +
                                (2.0 * t - 1.0) * end_point_m);
    }
  } else {
    const ScalarType denominator =
        (1.0 - t) * (1.0 - t) + 2.0 * weight_m * t * (1.0 - t) + t * t;
    PtBase<ScalarType> derivative =
        ((1.0 - t) * (1.0 - t) * start_point_m +
         2.0 * weight_m * t * (1.0 - t) * control_point_m +
         t * t * end_point_m);
    return derivative / denominator;
  }
}

template <class ScalarType>
inline PtBase<ScalarType> RationalBezierArcBase<ScalarType>::derivative(
    const ScalarType t) const {
  /* Defining constants and types */
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType FOUR = static_cast<ScalarType>(4);

  assert(t >= ZERO && t <= ONE);
  if (weight_m > static_cast<ScalarType>(1.0e15)) {
    if (t < ONE / TWO) {
      return PtBase<ScalarType>(TWO * (control_point_m - start_point_m));
    } else {
      return PtBase<ScalarType>(TWO * (end_point_m - control_point_m));
    }
  } else {
    ScalarType denominator =
        (ONE - t) * (ONE - t) + TWO * weight_m * t * (ONE - t) + t * t;
    denominator *= denominator;
    PtBase<ScalarType> derivative =
        (TWO * (start_point_m - end_point_m) * (ONE - weight_m) * t * t +
         FOUR * (start_point_m - control_point_m) * weight_m * t -
         TWO * (start_point_m - end_point_m) * t -
         TWO * weight_m * (start_point_m - control_point_m));
    return derivative / denominator;
  }
}

template <class ScalarType>
inline ScalarType RationalBezierArcBase<ScalarType>::arc_length(void) const {
  /* Defining constants and types */
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType THREE = static_cast<ScalarType>(3);
  const ScalarType FIVE = static_cast<ScalarType>(5);
  const ScalarType EIGHT = static_cast<ScalarType>(8);
  const ScalarType EIGHTEEN = static_cast<ScalarType>(18);

  // 3-point quadrature rule for arc_length calculation
  const ScalarType t0 = (ONE - sqrt(THREE / FIVE)) / TWO;
  const ScalarType t1 = ONE / TWO;
  const ScalarType t2 = (ONE + sqrt(THREE / FIVE)) / TWO;
  const ScalarType w0 = FIVE / EIGHTEEN;
  const ScalarType w1 = EIGHT / EIGHTEEN;
  const ScalarType w2 = FIVE / EIGHTEEN;
  PtBase<ScalarType> pt_0 = derivative(t0);
  PtBase<ScalarType> pt_1 = derivative(t1);
  PtBase<ScalarType> pt_2 = derivative(t2);
  const ScalarType norm0 =
      sqrt(pt_0[0] * pt_0[0] + pt_0[1] * pt_0[1] + pt_0[2] * pt_0[2]);
  const ScalarType norm1 =
      sqrt(pt_1[0] * pt_1[0] + pt_1[1] * pt_1[1] + pt_1[2] * pt_1[2]);
  const ScalarType norm2 =
      sqrt(pt_2[0] * pt_2[0] + pt_2[1] * pt_2[1] + pt_2[2] * pt_2[2]);
  return w0 * norm0 + w1 * norm1 + w2 * norm2;
}

template <class ScalarType>
inline std::uintptr_t RationalBezierArcBase<ScalarType>::start_point_id(
    void) const {
  return start_point_id_m;
}

template <class ScalarType>
inline std::uintptr_t RationalBezierArcBase<ScalarType>::end_point_id(
    void) const {
  return end_point_id_m;
}

template <class ScalarType>
inline RationalBezierArcBase<ScalarType>&
RationalBezierArcBase<ScalarType>::operator+=(const PtBase<ScalarType>& a_rhs) {
  this->start_point() += a_rhs;
  this->control_point() += a_rhs;
  this->end_point() += a_rhs;
  return (*this);
}

template <class ScalarType>
inline RationalBezierArcBase<ScalarType>
RationalBezierArcBase<ScalarType>::operator-(void) const {
  return RationalBezierArcBase(end_point_m, control_point_m, start_point_m,
                               weight_m);
}

template <class ScalarType>
inline ScalarType& RationalBezierArcBase<ScalarType>::weight(void) {
  return weight_m;
}

template <class ScalarType>
inline PtBase<ScalarType>& RationalBezierArcBase<ScalarType>::start_point(
    void) {
  return start_point_m;
}

template <class ScalarType>
inline PtBase<ScalarType>& RationalBezierArcBase<ScalarType>::control_point(
    void) {
  return control_point_m;
}

template <class ScalarType>
inline PtBase<ScalarType>& RationalBezierArcBase<ScalarType>::end_point(void) {
  return end_point_m;
}

template <class ScalarType>
inline std::uintptr_t RationalBezierArcBase<ScalarType>::start_point_id(void) {
  return start_point_id_m;
}

template <class ScalarType>
inline std::uintptr_t RationalBezierArcBase<ScalarType>::end_point_id(void) {
  return end_point_id_m;
}

template <class ScalarType>
inline void RationalBezierArcBase<ScalarType>::reset_start_point_id(
    std::uintptr_t a_id) {
  start_point_id_m = a_id;
}

template <class ScalarType>
inline void RationalBezierArcBase<ScalarType>::reset_end_point_id(
    std::uintptr_t a_id) {
  end_point_id_m = a_id;
}

template <class PtTypeWithGradient, class ScalarType>
inline RationalBezierArcWithGradientBase<
    PtTypeWithGradient, ScalarType>::RationalBezierArcWithGradientBase(void) {
  start_point_m = PtTypeWithGradient();
  control_point_m = PtTypeWithGradient();
  end_point_m = PtTypeWithGradient();
  weight_m = static_cast<ScalarType>(0);
  weight_gradient_m = PtTypeWithGradient::gradient_type();
}

template <class PtTypeWithGradient, class ScalarType>
inline RationalBezierArcWithGradientBase<PtTypeWithGradient, ScalarType>::
    RationalBezierArcWithGradientBase(
        const PtTypeWithGradient& a_start_pt,
        const PtTypeWithGradient& a_start_tangent,
        const PtTypeWithGradient& a_end_pt,
        const PtTypeWithGradient& a_end_tangent,
        const PtTypeWithGradient& a_plane_normal,
        const AlignedParaboloidBase<ScalarType>& a_paraboloid) {
  /* Defining constants and types */
  using gradient_type = typename PtTypeWithGradient::gradient_type;
  const ScalarType EPSILON = machine_epsilon<ScalarType>();
  const ScalarType ZERO = static_cast<ScalarType>(0);
  const ScalarType ONE = static_cast<ScalarType>(1);
  const ScalarType TWO = static_cast<ScalarType>(2);
  const ScalarType HALF = ONE / TWO;
  const ScalarType ONEHUNDRED = static_cast<ScalarType>(100);

  /* Function */
  start_point_m = a_start_pt;
  end_point_m = a_end_pt;
  const PtBase<ScalarType>& start_point = a_start_pt.getPt();
  const auto& start_point_grad = a_start_pt.getData();
  const PtBase<ScalarType>& start_tangent = a_start_tangent.getPt();
  const auto& start_tangent_grad = a_start_tangent.getData();
  const PtBase<ScalarType>& end_point = a_end_pt.getPt();
  const auto& end_point_grad = a_end_pt.getData();
  const PtBase<ScalarType>& end_tangent = a_end_tangent.getPt();
  const auto& end_tangent_grad = a_end_tangent.getData();

  // Compute control point
  const auto edge_vector_withgrad = end_point_m - start_point_m;
  const Normal edge_vector = Normal::fromPt(edge_vector_withgrad.getPt());
  const auto& edge_vector_grad = edge_vector_withgrad.getData();
  const auto average_pt_withgrad = HALF * (start_point_m + end_point_m);
  const auto& average_pt = average_pt_withgrad.getPt();
  const auto& average_pt_grad = average_pt_withgrad.getData();
  const auto& plane_normal = a_plane_normal.getPt();
  const auto& plane_normal_grad = a_plane_normal.getData();
  const PtBase<ScalarType> n_cross_t0 = PtBase<ScalarType>(
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
  // assert(std::fabs(Normal::fromPt(n_cross_t0) *
  // Normal::fromPt(end_tangent))
  // >
  //        EPSILON);
  const ScalarType n_cross_t0_dot_edge_vector = n_cross_t0[0] * edge_vector[0] +
                                                n_cross_t0[1] * edge_vector[1] +
                                                n_cross_t0[2] * edge_vector[2];
  const gradient_type n_cross_t0_dot_edge_vector_grad =
      n_cross_t0_grad[0] * edge_vector[0] +
      n_cross_t0[0] * edge_vector_grad[0] +
      n_cross_t0_grad[1] * edge_vector[1] +
      n_cross_t0[1] * edge_vector_grad[1] +
      n_cross_t0_grad[2] * edge_vector[2] + n_cross_t0[2] * edge_vector_grad[2];
  const ScalarType n_cross_t0_dot_end_tangent = n_cross_t0[0] * end_tangent[0] +
                                                n_cross_t0[1] * end_tangent[1] +
                                                n_cross_t0[2] * end_tangent[2];
  const gradient_type n_cross_t0_dot_end_tangent_grad =
      n_cross_t0_grad[0] * end_tangent[0] +
      n_cross_t0[0] * end_tangent_grad[0] +
      n_cross_t0_grad[1] * end_tangent[1] +
      n_cross_t0[1] * end_tangent_grad[1] +
      n_cross_t0_grad[2] * end_tangent[2] + n_cross_t0[2] * end_tangent_grad[2];
  const ScalarType lambda_1 =
      -(n_cross_t0_dot_edge_vector) / safelyEpsilon(n_cross_t0_dot_end_tangent);
  const gradient_type lambda_1_grad =
      -(n_cross_t0_dot_edge_vector_grad * n_cross_t0_dot_end_tangent -
        n_cross_t0_dot_edge_vector * n_cross_t0_dot_end_tangent_grad) /
      safelyEpsilon(n_cross_t0_dot_end_tangent * n_cross_t0_dot_end_tangent);
  const PtBase<ScalarType> control_point = end_point + lambda_1 * end_tangent;
  control_point_m = PtTypeWithGradient(control_point);
  auto& control_point_grad = control_point_m.getData();
  for (UnsignedIndex_t d = 0; d < 3; ++d) {
    control_point_grad[d] = end_point_grad[d] + lambda_1 * end_tangent_grad[d] +
                            lambda_1_grad * end_tangent[d];
  }
  const auto mid_to_control_withgrad = control_point_m - average_pt_withgrad;
  const auto projected_pt_withgrad =
      projectPtAlongHalfLineOntoParaboloidWithGradient<PtTypeWithGradient,
                                                       ScalarType>(
          a_paraboloid, mid_to_control_withgrad, average_pt_withgrad);
  const auto& projected_pt = projected_pt_withgrad.getPt();
  const auto& projected_pt_grad = projected_pt_withgrad.getData();
  weight_m = ONE;
  if (squaredMagnitude(projected_pt - control_point) < EPSILON * EPSILON) {
    weight_m = static_cast<ScalarType>(DBL_MAX);
    weight_gradient_m = gradient_type();
  } else {
    ScalarType previous_best = -static_cast<ScalarType>(DBL_MAX);
    UnsignedIndex_t chosen_id = 0;
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      const ScalarType denominator = projected_pt[d] - control_point[d];
      if (fabs(denominator) > previous_best) {
        chosen_id = d;
        previous_best = fabs(denominator);
        weight_m = (start_point[d] + end_point[d] - TWO * projected_pt[d]) /
                   safelyEpsilon(TWO * denominator);
      }
    }
    const ScalarType denominator =
        projected_pt[chosen_id] - control_point[chosen_id];
    const gradient_type den_grad =
        projected_pt_grad[chosen_id] - control_point_grad[chosen_id];
    weight_gradient_m =
        ((start_point_grad[chosen_id] + end_point_grad[chosen_id] -
          TWO * projected_pt_grad[chosen_id]) *
             denominator -
         (start_point[chosen_id] + end_point[chosen_id] -
          TWO * projected_pt[chosen_id]) *
             den_grad) /
        safelyEpsilon(TWO * denominator * denominator);
  }
}

template <class PtTypeWithGradient, class ScalarType>
inline RationalBezierArcBase<ScalarType>
RationalBezierArcWithGradientBase<PtTypeWithGradient, ScalarType>::arc(
    void) const {
  return RationalBezierArcBase(start_point_m.getPt(), control_point_m.getPt(),
                               end_point_m.getPt(), weight_m);
}

template <class PtTypeWithGradient, class ScalarType>
inline ScalarType
RationalBezierArcWithGradientBase<PtTypeWithGradient, ScalarType>::weight(
    void) const {
  return weight_m;
}

template <class PtTypeWithGradient, class ScalarType>
inline const typename PtTypeWithGradient::gradient_type&
RationalBezierArcWithGradientBase<PtTypeWithGradient,
                                  ScalarType>::weight_gradient(void) const {
  return weight_gradient_m;
}

template <class PtTypeWithGradient, class ScalarType>
inline const PtTypeWithGradient&
RationalBezierArcWithGradientBase<PtTypeWithGradient, ScalarType>::start_point(
    void) const {
  return start_point_m;
}

template <class PtTypeWithGradient, class ScalarType>
inline const PtTypeWithGradient& RationalBezierArcWithGradientBase<
    PtTypeWithGradient, ScalarType>::control_point(void) const {
  return control_point_m;
}

template <class PtTypeWithGradient, class ScalarType>
inline const PtTypeWithGradient&
RationalBezierArcWithGradientBase<PtTypeWithGradient, ScalarType>::end_point(
    void) const {
  return end_point_m;
}

template <class PtTypeWithGradient, class ScalarType>
inline ScalarType& RationalBezierArcWithGradientBase<PtTypeWithGradient,
                                                     ScalarType>::weight(void) {
  return weight_m;
}

template <class PtTypeWithGradient, class ScalarType>
inline typename PtTypeWithGradient::gradient_type&
RationalBezierArcWithGradientBase<PtTypeWithGradient,
                                  ScalarType>::weight_gradient(void) {
  return weight_gradient_m;
}

template <class PtTypeWithGradient, class ScalarType>
inline PtTypeWithGradient&
RationalBezierArcWithGradientBase<PtTypeWithGradient, ScalarType>::start_point(
    void) {
  return start_point_m;
}

template <class PtTypeWithGradient, class ScalarType>
inline PtTypeWithGradient& RationalBezierArcWithGradientBase<
    PtTypeWithGradient, ScalarType>::control_point(void) {
  return control_point_m;
}

template <class PtTypeWithGradient, class ScalarType>
inline PtTypeWithGradient&
RationalBezierArcWithGradientBase<PtTypeWithGradient, ScalarType>::end_point(
    void) {
  return end_point_m;
}

template <class ScalarType>
inline std::ostream& operator<<(
    std::ostream& out,
    const RationalBezierArcBase<ScalarType>& a_rational_bezier_arc) {
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
