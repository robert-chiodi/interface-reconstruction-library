// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_PARAMETRIZED_SURFACE_H_
#define IRL_PARABOLOID_RECONSTRUCTION_PARAMETRIZED_SURFACE_H_

#include <vector>

#ifdef IRL_USE_TRIANGLE
#include "external/triangle/triangle.h"
#endif

#include "irl/paraboloid_reconstruction/aligned_paraboloid.h"
#include "irl/paraboloid_reconstruction/ellipse.h"
#include "irl/paraboloid_reconstruction/rational_bezier_arc.h"
#include "irl/surface_mesher/triangulated_surface.h"

namespace IRL {

class NoSurfaceOutput {
 public:
  NoSurfaceOutput(void) = default;

 private:
};

/// \brief Parametrized surface defined by coeffs A,B of paraboloid + list of
/// rational Bézier arcs
class ParametrizedSurfaceOutput {
 public:
  /// \brief Default constructor.
  ParametrizedSurfaceOutput(void) = default;
  ParametrizedSurfaceOutput(const AlignedParaboloid& a_paraboloid);
  ~ParametrizedSurfaceOutput(void) = default;

  RationalBezierArc& operator[](const UnsignedIndex_t a_index);
  const RationalBezierArc& operator[](const UnsignedIndex_t a_index) const;
  const std::vector<RationalBezierArc>::size_type size(void) const;

  const AlignedParaboloid& getParaboloid(void) const;
  std::vector<RationalBezierArc>& getArcs(void);
  std::vector<Pt>& getPts(void);
  void addArc(const RationalBezierArc& a_rational_bezier_arc);
  void addPt(const Pt& a_pt);
  void clearArcs(void);
  void clearPts(void);
  void clear(void);

  TriangulatedSurfaceOutput triangulate(
      const double a_length_scale = -1.0,
      const UnsignedIndex_t a_nsplit = 5) const;

 private:
  AlignedParaboloid paraboloid_m;
  std::vector<Pt> pt_from_bezier_split_m;
  std::vector<RationalBezierArc> arc_list_m;
};

inline std::ostream& operator<<(
    std::ostream& out, const ParametrizedSurfaceOutput& a_parametrized_surface);

}  // namespace IRL

#include "irl/paraboloid_reconstruction/parametrized_surface.tpp"

#endif  // IRL_PARABOLOID_RECONSTRUCTION_PARAMETRIZED_SURFACE_H_
