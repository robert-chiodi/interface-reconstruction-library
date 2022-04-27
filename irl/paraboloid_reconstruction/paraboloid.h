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
// extern "C" {
// #include <quadmath.h>
// }
#include <cassert>
#include <ostream>

#include "irl/data_structures/stack_vector.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/pt_with_data.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/graphs/un_directed_graph_node.h"
#include "irl/paraboloid_reconstruction/aligned_paraboloid.h"
#include "irl/paraboloid_reconstruction/gradient_paraboloid.h"
#include "irl/planar_reconstruction/joined_reconstructions.h"
#include "irl/planar_reconstruction/planar_localizer.h"
#include "irl/planar_reconstruction/reconstruction_link.h"

namespace IRL {

// Paraboloid represented by a datum,
// reference frame, and an AlignedParaboloid defined as
// z + c + a*x^2 + b*y^2 = 0
class Paraboloid {
 public:
  Paraboloid(void) = default;

  /// a_reference_frame should have the normal vector in a_reference_frame[2],
  /// and the tangent vectors corresponding to a_coef_a and a_coef_b in
  /// element 0 and 1, respectively.
  Paraboloid(const Pt& a_datum, const ReferenceFrame& a_reference_frame,
             const double a_coef_a, const double a_coef_b);

  static Paraboloid createAlwaysAbove(void);

  static Paraboloid createAlwaysBelow(void);

  void setDatum(const Pt& a_datum);
  void setReferenceFrame(const ReferenceFrame& a_reference_frame);
  void setAlignedParaboloid(const AlignedParaboloid& a_aligned_paraboloid);

  const Pt& getDatum(void) const;
  const ReferenceFrame& getReferenceFrame(void) const;
  const AlignedParaboloid& getAlignedParaboloid(void) const;

  /// Indicates that the intersection should actually be performed.
  void markAsRealReconstruction(void);

  /// Marks paraboloid as being above any polyhedron (so any polyhedron will
  /// be unclipped).
  void markAsAlwaysAbove(void);

  /// Marks paraboloid as being below any polyhedron (so any polyhedron will
  /// be clipped).
  void markAsAlwaysBelow(void);

  /// Whether the paraboloid has been set to be above any polyhedron.
  bool isAlwaysAbove(void) const;

  /// Whether the paraboloid has been set to be below any polyhedron.
  bool isAlwaysBelow(void) const;

  /// Paraboloid cannot be a flipped reconstruction. Add this for ease of use
  /// with other routines that usually take planar reconstructions.
  static constexpr bool isFlipped(void) { return false; }

  /// \brief Since localizers are always convex, never flip.
  static constexpr double flip(void) { return 1.0; }

  /// \brief Return if cutting for gas phase is needed.
  static constexpr bool isNotFlipped(void) { return true; }

  ~Paraboloid(void) = default;

 private:
  Pt datum_m;
  ReferenceFrame frame_m;
  AlignedParaboloid paraboloid_m;
  std::array<bool, 2> place_infinite_shortcut_m;
};

using LocalizedParaboloid = JoinedReconstructions<PlanarLocalizer, Paraboloid>;
using LocalizedParaboloidLink =
    ReconstructionLink<LocalizedParaboloid, UnDirectedGraphNode>;

inline Pt conicCenter(const Plane& a_plane,
                      const AlignedParaboloid& a_paraboloid);

inline Normal getParaboloidSurfaceNormal(const AlignedParaboloid& a_paraboloid,
                                         const Pt& a_pt);

template <class PtTypeWithGradient>
inline PtTypeWithGradient getParaboloidSurfaceNormalWithGradient(
    const AlignedParaboloid& a_paraboloid, const PtTypeWithGradient& a_pt);

inline StackVector<double, 2> solveQuadratic(const double a, const double b,
                                             const double c);

template <class GradientType>
inline StackVector<std::pair<double, GradientType>, 2>
solveQuadraticWithGradient(const double a, const double b, const double c,
                           const GradientType& a_grad,
                           const GradientType& b_grad,
                           const GradientType& c_grad);

inline Pt projectPtAlongLineOntoParaboloid(
    const AlignedParaboloid& a_paraboloid, const Normal& a_line,
    const Pt& a_starting_pt);

inline Pt projectPtAlongHalfLineOntoParaboloid(
    const AlignedParaboloid& a_paraboloid, const Normal& a_line,
    const Pt& a_starting_pt);

template <class PtTypeWithGradient>
inline PtTypeWithGradient projectPtAlongHalfLineOntoParaboloidWithGradient(
    const AlignedParaboloid& a_paraboloid, const PtTypeWithGradient& a_line,
    const PtTypeWithGradient& a_starting_pt);

inline std::ostream& operator<<(std::ostream& out,
                                const Paraboloid& a_paraboloid);

}  // namespace IRL

#include "irl/paraboloid_reconstruction/paraboloid.tpp"

#endif  // IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_H_
