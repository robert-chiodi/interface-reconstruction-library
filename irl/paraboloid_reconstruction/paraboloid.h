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
template <class ScalarType>
class ParaboloidBase {
 public:
  using value_type = ScalarType;
  ParaboloidBase(void);

  /// a_reference_frame should have the normal vector in a_reference_frame[2],
  /// and the tangent vectors corresponding to a_coef_a and a_coef_b in
  /// element 0 and 1, respectively.
  ParaboloidBase(const PtBase<ScalarType>& a_datum,
                 const ReferenceFrameBase<ScalarType>& a_reference_frame,
                 const ScalarType a_coef_a, const ScalarType a_coef_b);

  static ParaboloidBase createAlwaysAbove(void);

  static ParaboloidBase createAlwaysBelow(void);

  void setDatum(const PtBase<ScalarType>& a_datum);
  void setReferenceFrame(
      const ReferenceFrameBase<ScalarType>& a_reference_frame);
  void setAlignedParaboloid(
      const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid);

  const PtBase<ScalarType>& getDatum(void) const;
  const ReferenceFrameBase<ScalarType>& getReferenceFrame(void) const;
  const AlignedParaboloidBase<ScalarType>& getAlignedParaboloid(void) const;

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
  static constexpr ScalarType flip(void) { return static_cast<ScalarType>(1); }

  /// \brief Return if cutting for gas phase is needed.
  static constexpr bool isNotFlipped(void) { return true; }

  ~ParaboloidBase(void) = default;

 private:
  PtBase<ScalarType> datum_m;
  ReferenceFrameBase<ScalarType> frame_m;
  AlignedParaboloidBase<ScalarType> paraboloid_m;
  std::array<bool, 2> place_infinite_shortcut_m;
};

using Paraboloid = ParaboloidBase<double>;

template <class ScalarType>
using LocalizedParaboloid =
    JoinedReconstructions<PlanarLocalizer, ParaboloidBase<ScalarType>>;
// using LocalizedParaboloid = LocalizedParaboloidBase<double>;

template <class ScalarType>
using LocalizedParaboloidLink =
    ReconstructionLink<LocalizedParaboloid<ScalarType>, UnDirectedGraphNode>;
// using LocalizedParaboloidLink = LocalizedParaboloidLinkBase<double>;

template <class ScalarType>
inline PtBase<ScalarType> conicCenter(
    const PlaneBase<ScalarType>& a_plane,
    const AlignedParaboloidBase<ScalarType>& a_paraboloid);

template <class ScalarType>
inline NormalBase<ScalarType> getParaboloidSurfaceNormal(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtBase<ScalarType>& a_pt);

template <class PtTypeWithGradient, class ScalarType>
inline PtTypeWithGradient getParaboloidSurfaceNormalWithGradient(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtTypeWithGradient& a_pt);

template <class ScalarType>
inline StackVector<ScalarType, 2> solveQuadratic(const ScalarType a,
                                                 const ScalarType b,
                                                 const ScalarType c);

template <class GradientType, class ScalarType>
inline StackVector<std::pair<ScalarType, GradientType>, 2>
solveQuadraticWithGradient(const ScalarType a, const ScalarType b,
                           const ScalarType c, const GradientType& a_grad,
                           const GradientType& b_grad,
                           const GradientType& c_grad);

template <class ScalarType>
inline PtBase<ScalarType> projectPtAlongLineOntoParaboloid(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const NormalBase<ScalarType>& a_line,
    const PtBase<ScalarType>& a_starting_pt);

template <class ScalarType>
inline PtBase<ScalarType> projectPtAlongHalfLineOntoParaboloid(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const NormalBase<ScalarType>& a_line,
    const PtBase<ScalarType>& a_starting_pt);

template <class PtTypeWithGradient, class ScalarType>
inline PtTypeWithGradient projectPtAlongHalfLineOntoParaboloidWithGradient(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtTypeWithGradient& a_line, const PtTypeWithGradient& a_starting_pt);

template <class ScalarType>
inline std::ostream& operator<<(std::ostream& out,
                                const ParaboloidBase<ScalarType>& a_paraboloid);

}  // namespace IRL

#include "irl/paraboloid_reconstruction/paraboloid.tpp"

#endif  // IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_H_
