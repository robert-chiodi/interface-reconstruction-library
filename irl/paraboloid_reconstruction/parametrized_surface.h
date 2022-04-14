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

#include "irl/geometry/general/normal.h"
#include "irl/paraboloid_reconstruction/ellipse.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/paraboloid_reconstruction/rational_bezier_arc.h"
#include "irl/surface_mesher/triangulated_surface.h"

namespace IRL {

template <class MomentType, class SurfaceType>
class AddSurfaceOutput {
 public:
  using moment_type = MomentType;
  using surface_type = SurfaceType;

  AddSurfaceOutput(void) = default;

  MomentType& getMoments(void);
  const MomentType& getMoments(void) const;

  SurfaceType& getSurface(void);
  const SurfaceType& getSurface(void) const;

 private:
  MomentType volume_moments_m;
  SurfaceType surface_m;
};

template <class C>
struct has_paraboloid_surface : std::false_type {};

template <class C>
struct has_paraboloid_surface<const C> : has_paraboloid_surface<C> {};

template <class MomentType, class SurfaceType>
struct has_paraboloid_surface<AddSurfaceOutput<MomentType, SurfaceType>>
    : std::true_type {};

class NoSurfaceOutput {
 public:
  NoSurfaceOutput(void) = default;
  ~NoSurfaceOutput(void) = default;

 private:
};

/// \brief Parametrized surface defined by coeffs A,B of paraboloid + list of
/// rational Bézier arcs
class ParametrizedSurfaceOutput {
 public:
  /// \brief Default constructor.
  ParametrizedSurfaceOutput(void) = default;
  ParametrizedSurfaceOutput(const Paraboloid& a_paraboloid);

  ParametrizedSurfaceOutput(const ParametrizedSurfaceOutput& a_rhs) = delete;
  ParametrizedSurfaceOutput(ParametrizedSurfaceOutput&& a_rhs);

  ParametrizedSurfaceOutput& operator=(const ParametrizedSurfaceOutput& a_rhs) =
      delete;
  ParametrizedSurfaceOutput& operator=(ParametrizedSurfaceOutput&& a_rhs);

  void setParaboloid(const Paraboloid& a_paraboloid);

  RationalBezierArc& operator[](const UnsignedIndex_t a_index);
  const RationalBezierArc& operator[](const UnsignedIndex_t a_index) const;
  const std::vector<RationalBezierArc>::size_type size(void) const;

  const Paraboloid& getParaboloid(void) const;
  std::vector<RationalBezierArc>& getArcs(void);
  std::vector<Pt*>& getPts(void);
  void addArc(const RationalBezierArc& a_rational_bezier_arc);
  void addPt(Pt* a_pt);
  void clearArcs(void);
  void clearPts(void);
  void clear(void);
  inline double getSurfaceArea(void);
  inline double getMeanCurvatureIntegral(void);
  inline double getGaussianCurvatureIntegral(void);
  inline Normal getAverageNormal(void);
  inline double getAverageMeanCurvature(void);
  inline double getAverageGaussianCurvature(void);

  TriangulatedSurfaceOutput triangulate(
      const double a_length_scale = -1.0,
      const UnsignedIndex_t a_nsplit = 5) const;

  ~ParametrizedSurfaceOutput(void);

 private:
  bool knows_surface_area_m;
  double surface_area_m;
  bool knows_avg_normal_m;
  Normal avg_normal_m;
  bool knows_int_mean_curv_m;
  double int_mean_curv_m;
  bool knows_int_gaussian_curv_m;
  double int_gaussian_curv_m;
  Paraboloid paraboloid_m;
  std::vector<Pt*> pt_from_bezier_split_m;
  std::vector<RationalBezierArc> arc_list_m;
};

inline std::ostream& operator<<(
    std::ostream& out, const ParametrizedSurfaceOutput& a_parametrized_surface);

}  // namespace IRL

#include "irl/paraboloid_reconstruction/parametrized_surface.tpp"

#endif  // IRL_PARABOLOID_RECONSTRUCTION_PARAMETRIZED_SURFACE_H_
