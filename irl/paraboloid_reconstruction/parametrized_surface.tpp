// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_PARAMETRIZED_SURFACE_TPP_
#define IRL_PARABOLOID_RECONSTRUCTION_PARAMETRIZED_SURFACE_TPP_

#include <fstream>
#include <iomanip>

#include "external/NumericalIntegration/NumericalIntegration.h"

namespace IRL {

template <class MomentType, class SurfaceType>
MomentType& AddSurfaceOutput<MomentType, SurfaceType>::getMoments(void) {
  return volume_moments_m;
}

template <class MomentType, class SurfaceType>
const MomentType& AddSurfaceOutput<MomentType, SurfaceType>::getMoments(
    void) const {
  return volume_moments_m;
}

template <class MomentType, class SurfaceType>
SurfaceType& AddSurfaceOutput<MomentType, SurfaceType>::getSurface(void) {
  return surface_m;
}

template <class MomentType, class SurfaceType>
const SurfaceType& AddSurfaceOutput<MomentType, SurfaceType>::getSurface(
    void) const {
  return surface_m;
}

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput()
    : knows_surface_area_m{false},
      knows_avg_normal_m{false},
      knows_int_mean_curv_m{false},
      knows_int_gaussian_curv_m{false},
      length_scale_m{-1.0} {}

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput(
    const Paraboloid& a_paraboloid)
    : paraboloid_m{a_paraboloid},
      knows_surface_area_m{false},
      knows_avg_normal_m{false},
      knows_int_mean_curv_m{false},
      knows_int_gaussian_curv_m{false},
      length_scale_m{-1.0} {}

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput(
    ParametrizedSurfaceOutput&& a_rhs)
    : paraboloid_m(a_rhs.paraboloid_m),
      arc_list_m(std::move(a_rhs.arc_list_m)),
      knows_surface_area_m(a_rhs.knows_surface_area_m),
      surface_area_m(a_rhs.surface_area_m),
      knows_avg_normal_m(a_rhs.knows_avg_normal_m),
      avg_normal_m(a_rhs.avg_normal_m),
      knows_int_mean_curv_m(a_rhs.knows_int_mean_curv_m),
      int_mean_curv_m(a_rhs.int_mean_curv_m),
      knows_int_gaussian_curv_m(a_rhs.knows_int_gaussian_curv_m),
      int_gaussian_curv_m(a_rhs.int_gaussian_curv_m),
      length_scale_m{a_rhs.length_scale_m} {}

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput(
    const ParametrizedSurfaceOutput& a_rhs)
    : paraboloid_m(a_rhs.paraboloid_m),
      arc_list_m(a_rhs.arc_list_m),
      knows_surface_area_m(a_rhs.knows_surface_area_m),
      surface_area_m(a_rhs.surface_area_m),
      knows_avg_normal_m(a_rhs.knows_avg_normal_m),
      avg_normal_m(a_rhs.avg_normal_m),
      knows_int_mean_curv_m(a_rhs.knows_int_mean_curv_m),
      int_mean_curv_m(a_rhs.int_mean_curv_m),
      knows_int_gaussian_curv_m(a_rhs.knows_int_gaussian_curv_m),
      int_gaussian_curv_m(a_rhs.int_gaussian_curv_m),
      length_scale_m{a_rhs.length_scale_m} {}

inline ParametrizedSurfaceOutput& ParametrizedSurfaceOutput::operator=(
    ParametrizedSurfaceOutput&& a_rhs) {
  if (this != &a_rhs) {
    paraboloid_m = a_rhs.paraboloid_m;
    arc_list_m = std::move(a_rhs.arc_list_m);
    knows_surface_area_m = a_rhs.knows_surface_area_m;
    surface_area_m = a_rhs.surface_area_m;
    knows_avg_normal_m = a_rhs.knows_avg_normal_m;
    avg_normal_m = a_rhs.avg_normal_m;
    knows_int_mean_curv_m = a_rhs.knows_int_mean_curv_m;
    int_mean_curv_m = a_rhs.int_mean_curv_m;
    knows_int_gaussian_curv_m = a_rhs.knows_int_gaussian_curv_m;
    int_gaussian_curv_m = a_rhs.int_gaussian_curv_m;
    length_scale_m = a_rhs.length_scale_m;
  }
  return *this;
}

inline ParametrizedSurfaceOutput& ParametrizedSurfaceOutput::operator=(
    const ParametrizedSurfaceOutput& a_rhs) {
  if (this != &a_rhs) {
    paraboloid_m = a_rhs.paraboloid_m;
    arc_list_m = a_rhs.arc_list_m;
    knows_surface_area_m = a_rhs.knows_surface_area_m;
    surface_area_m = a_rhs.surface_area_m;
    knows_avg_normal_m = a_rhs.knows_avg_normal_m;
    avg_normal_m = a_rhs.avg_normal_m;
    knows_int_mean_curv_m = a_rhs.knows_int_mean_curv_m;
    int_mean_curv_m = a_rhs.int_mean_curv_m;
    knows_int_gaussian_curv_m = a_rhs.knows_int_gaussian_curv_m;
    int_gaussian_curv_m = a_rhs.int_gaussian_curv_m;
    length_scale_m = a_rhs.length_scale_m;
  }
  return *this;
}

inline void ParametrizedSurfaceOutput::setLengthScale(
    const double a_length_scale) {
  length_scale_m = a_length_scale;
}

inline void ParametrizedSurfaceOutput::setParaboloid(
    const Paraboloid& a_paraboloid) {
  paraboloid_m = a_paraboloid;
}

inline RationalBezierArc& ParametrizedSurfaceOutput::operator[](
    const UnsignedIndex_t a_index) {
  return arc_list_m[a_index];
}

inline const RationalBezierArc& ParametrizedSurfaceOutput::operator[](
    const UnsignedIndex_t a_index) const {
  return arc_list_m[a_index];
}

inline const Paraboloid& ParametrizedSurfaceOutput::getParaboloid(void) const {
  return paraboloid_m;
}

inline std::vector<RationalBezierArc>& ParametrizedSurfaceOutput::getArcs(
    void) {
  return arc_list_m;
}

inline std::vector<Pt*>& ParametrizedSurfaceOutput::getPts(void) {
  return pt_from_bezier_split_m;
}

inline void ParametrizedSurfaceOutput::addArc(
    const RationalBezierArc& a_rational_bezier_arc) {
  arc_list_m.push_back(a_rational_bezier_arc);
}

inline void ParametrizedSurfaceOutput::addPt(Pt* a_pt) {
  pt_from_bezier_split_m.push_back(a_pt);
}

inline const std::vector<RationalBezierArc>::size_type
ParametrizedSurfaceOutput::size(void) const {
  return arc_list_m.size();
}

inline void ParametrizedSurfaceOutput::clearArcs(void) { arc_list_m.clear(); }

inline void ParametrizedSurfaceOutput::clearPts(void) {
  for (auto& elem : pt_from_bezier_split_m) {
    delete elem;
  }
  pt_from_bezier_split_m.clear();
}

inline void ParametrizedSurfaceOutput::clear(void) {
  this->clearArcs();
  this->clearPts();
}

inline ParametrizedSurfaceOutput::~ParametrizedSurfaceOutput(void) {
  for (auto elem : pt_from_bezier_split_m) {
    delete elem;
  }
}

class ArcContributionToSurfaceArea_Functor {
 public:
  ArcContributionToSurfaceArea_Functor(const RationalBezierArc& a_arc,
                                       const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else if (std::fabs(a) > std::fabs(b)) {
        const double primitive0 =
            (2. * a * pt0[0] *
                 std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
                           4. * (b * b) * (pt0[1] * pt0[1])) +
             (1. + 4. * (b * b) * (pt0[1] * pt0[1])) *
                 std::log(a * (2. * a * pt0[0] +
                               std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
                                         4. * (b * b) * (pt0[1] * pt0[1]))))) /
            (4. * a);
        const double primitive1 =
            (2. * a * pt1[0] *
                 std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
                           4. * (b * b) * (pt1[1] * pt1[1])) +
             (1. + 4. * (b * b) * (pt1[1] * pt1[1])) *
                 std::log(a * (2. * a * pt1[0] +
                               std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
                                         4. * (b * b) * (pt1[1] * pt1[1]))))) /
            (4. * a);
        return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
      } else {
        const double primitive0 =
            -(2. * b * pt0[1] *
                  std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
                            4. * (b * b) * (pt0[1] * pt0[1])) +
              (1. + 4. * (a * a) * (pt0[0] * pt0[0])) *
                  std::log(b *
                           (2. * b * pt0[1] +
                            std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0]) +
                                      4. * (b * b) * (pt0[1] * pt0[1]))))) /
            (4. * b);
        const double primitive1 =
            -(2. * b * pt1[1] *
                  std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
                            4. * (b * b) * (pt1[1] * pt1[1])) +
              (1. + 4. * (a * a) * (pt1[0] * pt1[0])) *
                  std::log(b *
                           (2. * b * pt1[1] +
                            std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0]) +
                                      4. * (b * b) * (pt1[1] * pt1[1]))))) /
            (4. * b);
        return 0.5 * (primitive0 * der0[0] + primitive1 * der1[0]);
      }
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return pt[0] * der[1];
      } else if (std::fabs(a) > std::fabs(b)) {
        const double primitive =
            (2. * pt[0] *
                 std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0]) +
                           4. * (b * b) * (pt[1] * pt[1])) +
             (1. + 4. * (b * b) * (pt[1] * pt[1])) *
                 std::log(-2. * std::fabs(a) * pt[0] +
                          std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0]) +
                                    4. * (b * b) * (pt[1] * pt[1]))) /
                 std::fabs(a)) /
            4.;
        if (std::isnan(primitive)) {
          std::cout << "Pr = " << pt << std::endl;
          std::cout << "Der = " << der << std::endl;
          std::cout << "a = " << a << std::endl;
          std::cout << "b = " << b << std::endl;
          std::cout << "Arc: weight = " << arc_m.weight() << std::endl;
          std::cout << "Arc: start = " << arc_m.start_point() << std::endl;
          std::cout << "Arc: ctrl  = " << arc_m.control_point() << std::endl;
          std::cout << "Arc: end   = " << arc_m.end_point() << std::endl;
          std::cout << "Primitive is NaN" << std::endl;
          exit(1);
        }
        if (std::isnan(der[1])) {
          std::cout << "der[1] is NaN" << std::endl;
          exit(1);
        }

        return primitive * der[1];
      } else {
        const double primitive =
            -(2. * pt[1] *
                  std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0]) +
                            4. * (b * b) * (pt[1] * pt[1])) +
              (1. + 4. * (a * a) * (pt[0] * pt[0])) *
                  std::log(-2. * std::fabs(b) * pt[1] +
                           std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0]) +
                                     4. * (b * b) * (pt[1] * pt[1]))) /
                  std::fabs(b)) /
            (4.);
        if (std::isnan(primitive)) {
          std::cout << "Pr = " << pt << std::endl;
          std::cout << "Der = " << der << std::endl;
          std::cout << "a = " << a << std::endl;
          std::cout << "b = " << b << std::endl;
          std::cout << "Arc: weight = " << arc_m.weight() << std::endl;
          std::cout << "Arc: start = " << arc_m.start_point() << std::endl;
          std::cout << "Arc: ctrl  = " << arc_m.control_point() << std::endl;
          std::cout << "Arc: end   = " << arc_m.end_point() << std::endl;
          std::cout << "Primitive is NaN" << std::endl;
          exit(1);
        }
        if (std::isnan(der[0])) {
          std::cout << "der[0] is NaN" << std::endl;
          exit(1);
        }

        return primitive * der[0];
      }
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};

class ArcContributionToNormalX_Functor {
 public:
  ArcContributionToNormalX_Functor(const RationalBezierArc& a_arc,
                                   const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      const double primitive0 = a * pt0[0] * pt0[0];
      const double primitive1 = a * pt1[0] * pt1[0];
      return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      const double primitive = a * pt[0] * pt[0];
      return primitive * der[1];
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};

class ArcContributionToNormalY_Functor {
 public:
  ArcContributionToNormalY_Functor(const RationalBezierArc& a_arc,
                                   const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      const double primitive0 = -b * pt0[1] * pt0[1];
      const double primitive1 = -b * pt1[1] * pt1[1];
      return 0.5 * (primitive0 * der0[0] + primitive1 * der1[0]);
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      const double primitive = -b * pt[1] * pt[1];
      return primitive * der[0];
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};

class ArcContributionToNormalZ_Functor {
 public:
  ArcContributionToNormalZ_Functor(const RationalBezierArc& a_arc,
                                   const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double primitive0 = pt0[0];
      const double primitive1 = pt1[0];
      return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double primitive = pt[0];
      return primitive * der[1];
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};

class ArcContributionToMeanCurvature_Functor {
 public:
  ArcContributionToMeanCurvature_Functor(const RationalBezierArc& a_arc,
                                         const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else if (std::fabs(a) > std::fabs(b)) {
        const double primitive0 =
            2. * b * pt0[0] +
            ((a + 4. * a * (b * b) * (pt0[1] * pt0[1]) -
              4. * (b * b * b) * (pt0[1] * pt0[1])) *
             std::atan((2. * a * pt0[0]) /
                       std::sqrt(1. + 4. * (b * b) * (pt0[1] * pt0[1])))) /
                (a * std::sqrt(1. + 4. * (b * b) * (pt0[1] * pt0[1])));
        const double primitive1 =
            2. * b * pt1[0] +
            ((a + 4. * a * (b * b) * (pt1[1] * pt1[1]) -
              4. * (b * b * b) * (pt1[1] * pt1[1])) *
             std::atan((2. * a * pt1[0]) /
                       std::sqrt(1. + 4. * (b * b) * (pt1[1] * pt1[1])))) /
                (a * std::sqrt(1. + 4. * (b * b) * (pt1[1] * pt1[1])));
        return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
      } else {
        const double primitive0 =
            -2. * a * pt0[1] -
            ((b - 4. * (a * a * a) * (pt0[0] * pt0[0]) +
              4. * (a * a) * b * (pt0[0] * pt0[0])) *
             std::atan((2. * b * pt0[1]) /
                       std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0])))) /
                (b * std::sqrt(1. + 4. * (a * a) * (pt0[0] * pt0[0])));
        const double primitive1 =
            -2. * a * pt1[1] -
            ((b - 4. * (a * a * a) * (pt1[0] * pt1[0]) +
              4. * (a * a) * b * (pt1[0] * pt1[0])) *
             std::atan((2. * b * pt1[1]) /
                       std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0])))) /
                (b * std::sqrt(1. + 4. * (a * a) * (pt1[0] * pt1[0])));
        return 0.5 * (primitive0 * der0[0] + primitive1 * der1[0]);
      }
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else if (std::fabs(a) > std::fabs(b)) {
        const double primitive =
            2. * b * pt[0] +
            ((a + 4. * a * (b * b) * (pt[1] * pt[1]) -
              4. * (b * b * b) * (pt[1] * pt[1])) *
             std::atan((2. * a * pt[0]) /
                       std::sqrt(1. + 4. * (b * b) * (pt[1] * pt[1])))) /
                (a * std::sqrt(1. + 4. * (b * b) * (pt[1] * pt[1])));
        return primitive * der[1];
      } else {
        const double primitive =
            -2. * a * pt[1] -
            ((b - 4. * (a * a * a) * (pt[0] * pt[0]) +
              4. * (a * a) * b * (pt[0] * pt[0])) *
             std::atan((2. * b * pt[1]) /
                       std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0])))) /
                (b * std::sqrt(1. + 4. * (a * a) * (pt[0] * pt[0])));
        return primitive * der[0];
      }
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};  // namespace IRL

class ArcContributionToGaussianCurvature_Functor {
 public:
  ArcContributionToGaussianCurvature_Functor(
      const RationalBezierArc& a_arc, const AlignedParaboloid& a_paraboloid)
      : arc_m(a_arc), paraboloid_m(a_paraboloid) {}

  double operator()(double a_t) const {
    const auto weight = arc_m.weight();
    if (weight > 1.0e15) {
      const auto pt0 = arc_m.point(0.5 * a_t);
      const auto pt1 = arc_m.point(0.5 + 0.5 * a_t);
      const auto der0 = arc_m.derivative(0.5 * a_t);
      const auto der1 = arc_m.derivative(0.5 + 0.5 * a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else {
        const double primitive0 =
            4.0 * a * b * pt0[0] /
            ((1.0 + 4.0 * b * b * pt0[1] * pt0[1]) *
             std::sqrt(1.0 + 4.0 * a * a * pt0[0] * pt0[0] +
                       4.0 * b * b * pt0[1] * pt0[1]));
        const double primitive1 =
            4.0 * a * b * pt1[0] /
            ((1.0 + 4.0 * b * b * pt1[1] * pt1[1]) *
             std::sqrt(1.0 + 4.0 * a * a * pt1[0] * pt1[0] +
                       4.0 * b * b * pt1[1] * pt1[1]));
        return 0.5 * (primitive0 * der0[1] + primitive1 * der1[1]);
      }
    } else {
      const auto pt = arc_m.point(a_t);
      const auto der = arc_m.derivative(a_t);
      const double a = paraboloid_m.a();
      const double b = paraboloid_m.b();
      if (std::fabs(a) < 10.0 * DBL_EPSILON &&
          std::fabs(b) < 10.0 * DBL_EPSILON) {
        return 0.0;
      } else {
        const double primitive = 4.0 * a * b * pt[0] /
                                 ((1.0 + 4.0 * b * b * pt[1] * pt[1]) *
                                  std::sqrt(1.0 + 4.0 * a * a * pt[0] * pt[0] +
                                            4.0 * b * b * pt[1] * pt[1]));
        return primitive * der[1];
      }
    }
  }

 private:
  const RationalBezierArc& arc_m;
  const AlignedParaboloid& paraboloid_m;
};  // namespace IRL

inline double ParametrizedSurfaceOutput::getSurfaceArea(void) {
  if (!knows_surface_area_m) {
    const UnsignedIndex_t nArcs = this->size();
    surface_area_m = 0.0;
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToSurfaceArea_Functor functor(arc_list_m[t],
                                                   aligned_paraboloid);

      // Define the integrator.
      Eigen::Integrator<double> integrator(limit);

      // Define a quadrature rule.
      Eigen::Integrator<double>::QuadratureRule quadrature_rule =
          Eigen::Integrator<double>::GaussKronrod61;

      // Integrate.
      surface_area_m += integrator.quadratureAdaptive(functor, 0.0, 1.0, epsabs,
                                                      epsrel, quadrature_rule);
    }
    knows_surface_area_m = true;
  }
  return surface_area_m;
}

inline Normal ParametrizedSurfaceOutput::getAverageNormal(void) {
  if (!knows_avg_normal_m) {
    const UnsignedIndex_t nArcs = this->size();
    avg_normal_m = Normal();
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToNormalX_Functor functorx(arc_list_m[t],
                                                aligned_paraboloid);
      ArcContributionToNormalY_Functor functory(arc_list_m[t],
                                                aligned_paraboloid);
      ArcContributionToNormalZ_Functor functorz(arc_list_m[t],
                                                aligned_paraboloid);

      // Define the integrator.
      Eigen::Integrator<double> integrator(limit);

      // Define a quadrature rule.
      Eigen::Integrator<double>::QuadratureRule quadrature_rule =
          Eigen::Integrator<double>::GaussKronrod61;

      // Integrate.
      avg_normal_m[0] += integrator.quadratureAdaptive(
          functorx, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
      avg_normal_m[1] += integrator.quadratureAdaptive(
          functory, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
      avg_normal_m[2] += integrator.quadratureAdaptive(
          functorz, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
    }
    avg_normal_m.normalize();
    knows_avg_normal_m = true;
  }
  return avg_normal_m;
}

inline double ParametrizedSurfaceOutput::getMeanCurvatureIntegral(void) {
  if (!knows_int_mean_curv_m) {
    const UnsignedIndex_t nArcs = this->size();
    int_mean_curv_m = 0.0;
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToMeanCurvature_Functor functor(arc_list_m[t],
                                                     aligned_paraboloid);

      // Define the integrator.
      Eigen::Integrator<double> integrator(limit);

      // Define a quadrature rule.
      Eigen::Integrator<double>::QuadratureRule quadrature_rule =
          Eigen::Integrator<double>::GaussKronrod61;

      // Integrate.
      int_mean_curv_m += integrator.quadratureAdaptive(
          functor, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
    }
    knows_int_mean_curv_m = true;
  }
  return int_mean_curv_m;
}

inline double ParametrizedSurfaceOutput::getAverageMeanCurvature(void) {
  return this->getMeanCurvatureIntegral() / safelyTiny(this->getSurfaceArea());
}

inline double ParametrizedSurfaceOutput::getGaussianCurvatureIntegral(void) {
  if (!knows_int_gaussian_curv_m) {
    const UnsignedIndex_t nArcs = this->size();
    int_gaussian_curv_m = 0.0;
    size_t limit = 128;

    const double epsabs = 10.0 * DBL_EPSILON;
    const double epsrel = 0.0;
    auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();
    for (std::size_t t = 0; t < nArcs; ++t) {
      // Define the functor
      ArcContributionToGaussianCurvature_Functor functor(arc_list_m[t],
                                                         aligned_paraboloid);

      // Define the integrator.
      Eigen::Integrator<double> integrator(limit);

      // Define a quadrature rule.
      Eigen::Integrator<double>::QuadratureRule quadrature_rule =
          Eigen::Integrator<double>::GaussKronrod61;

      // Integrate.
      int_gaussian_curv_m += integrator.quadratureAdaptive(
          functor, 0.0, 1.0, epsabs, epsrel, quadrature_rule);
    }
    knows_int_gaussian_curv_m = true;
  }
  return int_gaussian_curv_m;
}

inline double ParametrizedSurfaceOutput::getAverageGaussianCurvature(void) {
  return this->getGaussianCurvatureIntegral() /
         safelyTiny(this->getSurfaceArea());
}

inline Normal ParametrizedSurfaceOutput::getNormalAligned(const Pt a_pt) {
  auto& aligned_paraboloid = this->getParaboloid().getAlignedParaboloid();
  auto aligned_normal = getParaboloidSurfaceNormal(aligned_paraboloid, a_pt);
  aligned_normal.normalize();
  return aligned_normal;
}

inline Normal ParametrizedSurfaceOutput::getNormalNonAligned(const Pt a_pt) {
  const auto& datum = this->getParaboloid().getDatum();
  const auto& ref_frame = this->getParaboloid().getReferenceFrame();
  // assert(ref_frame.isOrthonormalBasis());
  const Pt original_pt = a_pt - datum;
  auto aligned_pt = a_pt;
  for (std::size_t n = 0; n < 3; ++n) {
    aligned_pt[n] = ref_frame[n] * original_pt;
  }
  auto aligned_normal = this->getNormalAligned(aligned_pt);
  auto normal = Normal();
  for (std::size_t d = 0; d < 3; ++d) {
    for (std::size_t n = 0; n < 3; ++n) {
      normal[n] += ref_frame[d][n] * aligned_normal[d];
    }
  }
  return normal;
}

inline double ParametrizedSurfaceOutput::getMeanCurvatureAligned(
    const Pt a_pt) {
  auto& aligned_paraboloid = this->getParaboloid().getAlignedParaboloid();
  return (2. * (aligned_paraboloid.a() + aligned_paraboloid.b() +
                4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
                    aligned_paraboloid.b() * (a_pt[0] * a_pt[0]) +
                4. * aligned_paraboloid.a() *
                    (aligned_paraboloid.b() * aligned_paraboloid.b()) *
                    (a_pt[1] * a_pt[1]))) /
         std::pow(1. +
                      4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
                          (a_pt[0] * a_pt[0]) +
                      4. * (aligned_paraboloid.b() * aligned_paraboloid.b()) *
                          (a_pt[1] * a_pt[1]),
                  1.5);
}

inline double ParametrizedSurfaceOutput::getMeanCurvatureNonAligned(
    const Pt a_pt) {
  const auto& datum = this->getParaboloid().getDatum();
  const auto& ref_frame = this->getParaboloid().getReferenceFrame();
  // assert(ref_frame.isOrthonormalBasis());
  const Pt original_pt = a_pt - datum;
  auto aligned_pt = a_pt;
  for (std::size_t n = 0; n < 3; ++n) {
    aligned_pt[n] = ref_frame[n] * original_pt;
  }
  return this->getMeanCurvatureAligned(aligned_pt);
}

inline double ParametrizedSurfaceOutput::getGaussianCurvatureAligned(
    const Pt a_pt) {
  auto& aligned_paraboloid = this->getParaboloid().getAlignedParaboloid();
  return 4. * aligned_paraboloid.a() * aligned_paraboloid.b() /
         ((1. +
           4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
               (a_pt[0] * a_pt[0]) +
           4. * (aligned_paraboloid.b() * aligned_paraboloid.b()) *
               (a_pt[1] * a_pt[1])) *
          (1. +
           4. * (aligned_paraboloid.a() * aligned_paraboloid.a()) *
               (a_pt[0] * a_pt[0]) +
           4. * (aligned_paraboloid.b() * aligned_paraboloid.b()) *
               (a_pt[1] * a_pt[1])));
}

inline double ParametrizedSurfaceOutput::getGaussianCurvatureNonAligned(
    const Pt a_pt) {
  const auto& datum = this->getParaboloid().getDatum();
  const auto& ref_frame = this->getParaboloid().getReferenceFrame();
  // assert(ref_frame.isOrthonormalBasis());
  const Pt original_pt = a_pt - datum;
  auto aligned_pt = a_pt;
  for (std::size_t n = 0; n < 3; ++n) {
    aligned_pt[n] = ref_frame[n] * original_pt;
  }
  return this->getGaussianCurvatureAligned(aligned_pt);
}

inline TriangulatedSurfaceOutput ParametrizedSurfaceOutput::triangulate(
    const double a_length_scale, const UnsignedIndex_t a_nsplit) const {
#ifdef IRL_NO_USE_TRIANGLE
  std::cout << "IRL not compiled with the Triangle library. Exiting..."
            << std::endl;
  std::exit(-1);
#endif

  const UnsignedIndex_t nArcs = this->size();
  double length_scale, length_scale_ref = length_scale_m;
  if (a_length_scale > 0.0) {
    length_scale_ref = a_length_scale;
  }
  const auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();

  std::vector<std::vector<RationalBezierArc>> list_of_closed_curves;
  std::vector<bool> visited(nArcs, false);

  // First, we need to order the arcs so as to form closed curves
  double min_arc_length = DBL_MAX;
  for (std::size_t t = 0; t < nArcs; ++t) {
    if (visited[t]) {
      continue;
    }
    visited[t] = true;
    // Start with next available arc
    list_of_closed_curves.push_back(
        std::vector<RationalBezierArc>({arc_list_m[t]}));
    const std::uintptr_t start_id = arc_list_m[t].start_point_id();
    std::uintptr_t end_id = arc_list_m[t].end_point_id();
    while (end_id != start_id) {
      for (std::size_t e = t + 1; e < nArcs; ++e) {
        if (arc_list_m[e].start_point_id() == end_id) {
          visited[e] = true;
          list_of_closed_curves.back().push_back(arc_list_m[e]);
          end_id = arc_list_m[e].end_point_id();
          break;
        }
      }
    }
  }

  // Second, we approximate the arc length of the arc, so as to know how
  // many times it needs to be split
  std::vector<REAL> input_points;
  std::vector<REAL> input_holes;
  std::vector<int> input_segments;
  const UnsignedIndex_t nCurves =
      static_cast<UnsignedIndex_t>(list_of_closed_curves.size());
  // Loop over curves
  UnsignedIndex_t start_points = 0;
  for (UnsignedIndex_t i = 0; i < nCurves; ++i) {
    const UnsignedIndex_t nLocalArcs = list_of_closed_curves[i].size();
    // Loop over arcs of curve
    UnsignedIndex_t added_points = 0;
    double signed_area = 0.0;
    for (UnsignedIndex_t j = 0; j < nLocalArcs; ++j) {
      // Compute approximate arc length
      const RationalBezierArc& arc = list_of_closed_curves[i][j];
      const auto& sp = arc.start_point();
      const auto& ep = arc.start_point();
      const double arc_length = arc.arc_length();
      signed_area += (sp[0] * ep[1] - ep[0] * sp[1]);

      // Split arc
      UnsignedIndex_t nSplit = a_nsplit <= 0 ? 1 : a_nsplit;
      if (length_scale_ref > 0.0) {
        nSplit = static_cast<UnsignedIndex_t>(arc_length / length_scale_ref);
        nSplit = nSplit < a_nsplit ? a_nsplit : nSplit;
      }
      const double step = 1.0 / static_cast<double>(nSplit);
      length_scale = std::min(length_scale, step * arc_length);
      if (length_scale_ref > 0.0) length_scale = length_scale_ref;
      added_points += nSplit;
      const auto start_ind = input_points.size();
      input_points.resize(start_ind + 2 * nSplit);
      auto loc = input_points.begin() + start_ind;
      for (UnsignedIndex_t k = 1; k <= nSplit; ++k) {
        const double t = static_cast<double>(k) * step;
        const auto pt = arc.point(t);
        *(loc++) = pt[0];
        *(loc++) = pt[1];
      }
    }

    if (signed_area < 0.0) {
      // Add hole
      const auto p1x = input_points[start_points];
      const auto p1y = input_points[start_points + 1];
      const auto p2x = input_points[start_points + 2];
      const auto p2y = input_points[start_points + 3];
      std::array<double, 2> hole_location{
          {0.5 * (p1x + p2x), 0.5 * (p1y + p2y)}};
      Normal shift_dir = Normal(p2y - p1y, p1x - p2x, 0.0);
      shift_dir.normalize();
      const auto start_ind = input_points.size();
      input_holes.resize(start_ind + 2);
      input_holes[start_ind] =
          hole_location[0] + (5.0 * DBL_EPSILON) * shift_dir[0];
      input_holes[start_ind + 1] =
          hole_location[1] + (5.0 * DBL_EPSILON) * shift_dir[1];
    }

    // Create segments
    const int seg_size = input_segments.size();
    input_segments.resize(seg_size + 2 * (added_points));
    auto seg_loc = input_segments.begin() + seg_size;
    *(seg_loc++) = start_points + added_points - 1;
    *(seg_loc++) = start_points;
    for (UnsignedIndex_t j = start_points; j < start_points + added_points - 1;
         ++j) {
      *(seg_loc++) = j;
      *(seg_loc++) = j + 1;
    }
    start_points += added_points;
  }

  TriangulatedSurfaceOutput returned_surface;

  // Below section is for Triangle library
  if (input_points.size() > 0) {
    // Calling triangulation library
    struct triangulateio in = {0}, out = {0};
    in.numberofpoints = input_points.size() / 2;
    in.pointlist = input_points.data();

    std::vector<int> pointmarkerlist(in.numberofpoints, 1);
    in.pointmarkerlist = pointmarkerlist.data();

    in.numberofsegments = input_segments.size() / 2;
    in.segmentlist = input_segments.data();
    std::vector<int> segmentmarkerlist(in.numberofsegments, 1);
    in.segmentmarkerlist = segmentmarkerlist.data();

    in.numberofholes = input_holes.size() / 2;
    if (in.numberofholes > 0) {
      in.holelist = input_holes.data();
    }

    char flags[50];
    sprintf(flags, "pzqYYia%.15feQ", 0.5 * length_scale * length_scale);
    triangulate_from_lib(flags, &in, &out, (struct triangulateio*)NULL);

    auto& vlist = returned_surface.getVertexList();
    vlist.resize(out.numberofpoints);
    for (UnsignedIndex_t i = 0; i < out.numberofpoints; ++i) {
      const double x = out.pointlist[2 * i + 0];
      const double y = out.pointlist[2 * i + 1];
      const double z =
          -aligned_paraboloid.a() * x * x - aligned_paraboloid.b() * y * y;
      vlist[i] = Pt(x, y, z);
    }

    // Translate and rotate triangulated surface vertices
    const auto& datum = paraboloid_m.getDatum();
    const auto& ref_frame = paraboloid_m.getReferenceFrame();
    for (auto& vertex : vlist) {
      const Pt base_pt = vertex;
      vertex = Pt(0.0, 0.0, 0.0);
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          vertex[n] += ref_frame[d][n] * base_pt[d];
        }
      }
      vertex += datum;
    }

    for (UnsignedIndex_t i = 0; i < out.numberofedges; ++i) {
      if (out.edgemarkerlist[i] == 1) {
        returned_surface.addBoundaryEdge(out.edgelist[2 * i],
                                         out.edgelist[2 * i + 1]);
      }
    }

    auto& tlist = returned_surface.getTriangleList();
    tlist.resize(out.numberoftriangles,
                 TriangulatedSurfaceOutput::TriangleStorage::value_type::
                     fromNoExistencePlane(vlist, {0, 0, 0}));
    for (UnsignedIndex_t i = 0; i < out.numberoftriangles; ++i) {
      tlist[i] = TriangulatedSurfaceOutput::TriangleStorage::value_type::
          fromNoExistencePlane(
              vlist,
              {static_cast<UnsignedIndex_t>(out.trianglelist[3 * i + 0]),
               static_cast<UnsignedIndex_t>(out.trianglelist[3 * i + 1]),
               static_cast<UnsignedIndex_t>(out.trianglelist[3 * i + 2])});
    }

    /* free all allocated arrays, including those allocated by Triangle.
     */
    // free(in.pointlist);
    free(in.pointattributelist);
    // free(in.pointmarkerlist);
    free(in.trianglelist);
    free(in.triangleattributelist);
    free(in.trianglearealist);
    free(in.neighborlist);
    // free(in.segmentlist);
    // free(in.segmentmarkerlist);
    // free(in.holelist);
    free(in.regionlist);
    free(in.edgelist);
    free(in.edgemarkerlist);
    free(in.normlist);
    free(out.pointlist);
    free(out.pointattributelist);
    free(out.pointmarkerlist);
    free(out.trianglelist);
    free(out.triangleattributelist);
    free(out.trianglearealist);
    free(out.neighborlist);
    free(out.segmentlist);
    free(out.segmentmarkerlist);
    free(out.regionlist);
    free(out.edgelist);
    free(out.edgemarkerlist);
    free(out.normlist);
  }

  return returned_surface;
}

inline std::ostream& operator<<(
    std::ostream& out,
    const ParametrizedSurfaceOutput& a_parametrized_surface) {
  const auto& aligned_paraboloid =
      a_parametrized_surface.getParaboloid().getAlignedParaboloid();
  out.precision(16);
  out << std::scientific << aligned_paraboloid.a() << " "
      << aligned_paraboloid.b() << std::endl;
  for (UnsignedIndex_t i = 0; i < a_parametrized_surface.size(); ++i) {
    out << a_parametrized_surface[i];
    if (i < a_parametrized_surface.size() - 1) out << std::endl;
  }
  return out;
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTION_PARAMETRIZED_SURFACE_TPP_
