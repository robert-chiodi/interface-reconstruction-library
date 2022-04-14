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
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/graphs/un_directed_graph_node.h"
#include "irl/paraboloid_reconstruction/aligned_paraboloid.h"
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
  /// and the tangent vectors corresponding to a_coef_a and a_coef_b in element
  /// 0 and 1, respectively.
  Paraboloid(const Pt& a_datum, const ReferenceFrame& a_reference_frame,
             const double a_coef_a, const double a_coef_b);

  void setDatum(const Pt& a_datum);
  void setReferenceFrame(const ReferenceFrame& a_reference_frame);
  void setAlignedParaboloid(const AlignedParaboloid& a_aligned_paraboloid);

  const Pt& getDatum(void) const;
  const ReferenceFrame& getReferenceFrame(void) const;
  const AlignedParaboloid& getAlignedParaboloid(void) const;

  ~Paraboloid(void) = default;

 private:
  Pt datum_m;
  ReferenceFrame frame_m;
  AlignedParaboloid paraboloid_m;
};

using LocalizedParaboloid = JoinedReconstructions<PlanarLocalizer, Paraboloid>;
using LocalizedParaboloidLink =
    ReconstructionLink<LocalizedParaboloid, UnDirectedGraphNode>;

inline Normal getParaboloidSurfaceNormal(const AlignedParaboloid& a_paraboloid,
                                         const Pt& a_pt);

// Returns solution to quadratic equation solve.
// The smallest solution will always be first.
inline StackVector<double, 2> solveQuadratic(const double a, const double b,
                                             const double c);
inline Pt projectPtAlongLineOntoParaboloid(
    const AlignedParaboloid& a_paraboloid, const Normal& a_line,
    const Pt& a_starting_pt);

}  // namespace IRL

#include "irl/paraboloid_reconstruction/paraboloid.tpp"

#endif  // IRL_PARABOLOID_RECONSTRUCTIONS_PARABOLOID_H_
