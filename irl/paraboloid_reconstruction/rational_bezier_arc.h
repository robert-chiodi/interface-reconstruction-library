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
template <class ScalarType>
class RationalBezierArcBase {
 public:
  using value_type = ScalarType;
  /// \brief Default constructor.
  RationalBezierArcBase(void);

  /// \brief Constructor that initializes the rational Bèzier arc
  RationalBezierArcBase(const PtBase<ScalarType>& a_start_pt,
                        const PtBase<ScalarType>& a_control_pt,
                        const PtBase<ScalarType>& a_end_pt,
                        const ScalarType a_weight);

  /// \brief Constructor that initializes the rational Bèzier arc
  RationalBezierArcBase(const PtBase<ScalarType>& a_start_pt,
                        const PtBase<ScalarType>& a_control_pt,
                        const PtBase<ScalarType>& a_end_pt,
                        std::uintptr_t a_start_id, std::uintptr_t a_end_id,
                        const ScalarType a_weight);

  /// \brief Constructor that initializes the rational Bèzier arc by computing
  /// the rational weight
  RationalBezierArcBase(const PtBase<ScalarType>& a_start_pt,
                        const NormalBase<ScalarType>& a_start_tangent,
                        const PtBase<ScalarType>& a_end_pt,
                        const NormalBase<ScalarType>& a_end_tangent,
                        const NormalBase<ScalarType>& a_plane_normal,
                        const AlignedParaboloidBase<ScalarType>& a_paraboloid);

  /// \brief Return const weight.
  const ScalarType& weight(void) const;
  /// \brief Return const reference to stored start point.
  const PtBase<ScalarType>& start_point(void) const;
  /// \brief Return const reference to stored control point.
  const PtBase<ScalarType>& control_point(void) const;
  /// \brief Return const reference to stored end point.
  const PtBase<ScalarType>& end_point(void) const;
  /// \brief Return point evaluated at parameter a_t.
  PtBase<ScalarType> point(const ScalarType a_t) const;
  /// \brief Return derivative of the curve with respect to a_t.
  PtBase<ScalarType> derivative(const ScalarType a_t) const;
  /// \brief Return const approximation of arc_length.
  ScalarType arc_length(void) const;
  /// \brief Return const reference to stored start point address.
  std::uintptr_t start_point_id(void) const;
  /// \brief Return const reference to stored end point address.
  std::uintptr_t end_point_id(void) const;
  /// \brief Resets stored start point address.
  void reset_start_point_id(std::uintptr_t a_id);
  /// \brief Resets stored end point address.
  void reset_end_point_id(std::uintptr_t a_id);

  /// \brief Return const weight.
  ScalarType& weight(void);
  /// \brief Return const reference to stored start point.
  PtBase<ScalarType>& start_point(void);
  /// \brief Return const reference to stored control point.
  PtBase<ScalarType>& control_point(void);
  /// \brief Return const reference to stored end point.
  PtBase<ScalarType>& end_point(void);

  /// \brief Overload += operator to update moments.
  RationalBezierArcBase& operator+=(const PtBase<ScalarType>& a_rhs);
  // \brief Unary minus operator
  RationalBezierArcBase operator-(void) const;

  /// \brief Default destructor.
  ~RationalBezierArcBase(void) = default;

 private:
  /// \brief Return reference to stored start point address.
  std::uintptr_t start_point_id(void);
  /// \brief Return reference to stored end point address.
  std::uintptr_t end_point_id(void);

  PtBase<ScalarType> start_point_m;    // Start point.
  PtBase<ScalarType> control_point_m;  // Control point.
  PtBase<ScalarType> end_point_m;      // End point.
  std::uintptr_t start_point_id_m;     // Start point address.
  std::uintptr_t end_point_id_m;       // End point address.
  ScalarType weight_m;                 // Weight
};

template <class PtTypeWithGradient, class ScalarType>
/// \brief Rational Bézier arc defined by end points + control point +
/// weight
class RationalBezierArcWithGradientBase {
 public:
  using gradient_type = typename PtTypeWithGradient::gradient_type;
  using value_type = ScalarType;

  /// \brief Default constructor.
  RationalBezierArcWithGradientBase(void);

  /// \brief Constructor that initializes the rational Bèzier arc by computing
  /// the rational weight
  RationalBezierArcWithGradientBase(
      const PtTypeWithGradient& a_start_pt,
      const PtTypeWithGradient& a_start_tangent,
      const PtTypeWithGradient& a_end_pt,
      const PtTypeWithGradient& a_end_tangent,
      const PtTypeWithGradient& a_plane_normal,
      const AlignedParaboloidBase<ScalarType>& a_paraboloid);

  /// \brief Return arc.
  RationalBezierArcBase<ScalarType> arc(void) const;
  /// \brief Return const weight.
  ScalarType weight(void) const;
  /// \brief Return const weight.
  const gradient_type& weight_gradient(void) const;
  /// \brief Return const reference to stored start point.
  const PtTypeWithGradient& start_point(void) const;
  /// \brief Return const reference to stored control point.
  const PtTypeWithGradient& control_point(void) const;
  /// \brief Return const reference to stored end point.
  const PtTypeWithGradient& end_point(void) const;

  /// \brief Default destructor.
  ~RationalBezierArcWithGradientBase(void) = default;

 private:
  /// \brief Return const weight.
  ScalarType& weight(void);
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
  ScalarType weight_m;                 // Weight
  gradient_type weight_gradient_m;     // Weight
};

template <class ScalarType>
inline std::ostream& operator<<(
    std::ostream& out,
    const RationalBezierArcBase<ScalarType>& a_rational_bezier_arc);

using RationalBezierArc = RationalBezierArcBase<double>;

template <class PtTypeWithGradient>
using RationalBezierArcWithGradient =
    RationalBezierArcWithGradientBase<PtTypeWithGradient, double>;

}  // namespace IRL

#include "irl/paraboloid_reconstruction/rational_bezier_arc.tpp"

#endif  // IRL_PARABOLOID_RECONSTRUCTION_RATIONAL_BEZIER_ARC_H_
