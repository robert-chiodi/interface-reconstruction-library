// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_RATIONAL_BEZIER_ARC_H_
#define IRL_PARABOLOID_RECONSTRUCTION_RATIONAL_BEZIER_ARC_H_

#include <cstdint>
#include <utility>

#include "irl/geometry/general/pt.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

/// \brief Rational Bézier arc defined by end points + control point + weight
class RationalBezierArc {
 public:
  /// \brief Default constructor.
  RationalBezierArc(void);

  /// \brief Constructor that initializes the rational Bèzier arc
  RationalBezierArc(const Pt& a_start_pt, const Pt& a_control_pt,
                    const Pt& a_end_pt, const double a_weight);

  /// \brief Constructor that initializes the rational Bèzier arc by computing
  /// the rational weight
  RationalBezierArc(const Pt& a_start_pt, const Normal& a_start_tangent,
                    const Pt& a_end_pt, const Normal& a_end_tangent,
                    const Normal& a_plane_normal,
                    const AlignedParaboloid& a_paraboloid);

  /// \brief Return const weight.
  const double& weight(void) const;
  /// \brief Return const reference to stored start point.
  const Pt& start_point(void) const;
  /// \brief Return const reference to stored control point.
  const Pt& control_point(void) const;
  /// \brief Return const reference to stored end point.
  const Pt& end_point(void) const;
  /// \brief Return point evaluated at parameter a_t.
  Pt point(const double a_t) const;
  /// \brief Return derivative of the curve with respect to a_t.
  Pt derivative(const double a_t) const;
  /// \brief Return const approximation of arc_length.
  double arc_length(void) const;
  /// \brief Return const reference to stored start point address.
  std::uintptr_t start_point_id(void) const;
  /// \brief Return const reference to stored end point address.
  std::uintptr_t end_point_id(void) const;

  /// \brief Return const weight.
  double& weight(void);
  /// \brief Return const reference to stored start point.
  Pt& start_point(void);
  /// \brief Return const reference to stored control point.
  Pt& control_point(void);
  /// \brief Return const reference to stored end point.
  Pt& end_point(void);

  /// \brief Overload += operator to update moments.
  RationalBezierArc& operator+=(const Pt& a_rhs);
  // \brief Unary minus operator
  RationalBezierArc operator-(void) const;

  /// \brief Default destructor.
  ~RationalBezierArc(void) = default;

 private:
  /// \brief Return reference to stored start point address.
  std::uintptr_t start_point_id(void);
  /// \brief Return reference to stored end point address.
  std::uintptr_t end_point_id(void);

  Pt start_point_m;                 // Start point.
  Pt control_point_m;               // Control point.
  Pt end_point_m;                   // End point.
  std::uintptr_t start_point_id_m;  // Start point address.
  std::uintptr_t end_point_id_m;    // End point address.
  double weight_m;                  // Weight
};

template <class PtTypeWithGradient>
/// \brief Rational Bézier arc defined by end points + control point +
/// weight
class RationalBezierArcWithGradient {
 public:
  using gradient_type = typename PtTypeWithGradient::gradient_type;

  /// \brief Default constructor.
  RationalBezierArcWithGradient(void);

  /// \brief Constructor that initializes the rational Bèzier arc by computing
  /// the rational weight
  RationalBezierArcWithGradient(const PtTypeWithGradient& a_start_pt,
                                const PtTypeWithGradient& a_start_tangent,
                                const PtTypeWithGradient& a_end_pt,
                                const PtTypeWithGradient& a_end_tangent,
                                const PtTypeWithGradient& a_plane_normal,
                                const AlignedParaboloid& a_paraboloid);

  /// \brief Return arc.
  RationalBezierArc arc(void) const;
  /// \brief Return const weight.
  double weight(void) const;
  /// \brief Return const weight.
  const gradient_type& weight_gradient(void) const;
  /// \brief Return const reference to stored start point.
  const PtTypeWithGradient& start_point(void) const;
  /// \brief Return const reference to stored control point.
  const PtTypeWithGradient& control_point(void) const;
  /// \brief Return const reference to stored end point.
  const PtTypeWithGradient& end_point(void) const;

  /// \brief Default destructor.
  ~RationalBezierArcWithGradient(void) = default;

 private:
  /// \brief Return const weight.
  double& weight(void);
  /// \brief Return const weight.
  gradient_type& weight_gradient(void);
  /// \brief Return const reference to stored start point.
  PtTypeWithGradient& start_point(void);
  /// \brief Return const reference to stored control point.
  PtTypeWithGradient& control_point(void);
  /// \brief Return const reference to stored end point.
  PtTypeWithGradient& end_point(void);

  PtTypeWithGradient start_point_m;    // Start point.
  PtTypeWithGradient control_point_m;  // Control point.
  PtTypeWithGradient end_point_m;      // End point.
  double weight_m;                     // Weight
  gradient_type weight_gradient_m;     // Weight
};

inline std::ostream& operator<<(std::ostream& out,
                                const RationalBezierArc& a_rational_bezier_arc);

}  // namespace IRL

#include "irl/paraboloid_reconstruction/rational_bezier_arc.tpp"

#endif  // IRL_PARABOLOID_RECONSTRUCTION_RATIONAL_BEZIER_ARC_H_
