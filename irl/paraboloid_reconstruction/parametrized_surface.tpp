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

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput(
    const Paraboloid& a_paraboloid)
    : paraboloid_m{a_paraboloid} {}

inline ParametrizedSurfaceOutput::ParametrizedSurfaceOutput(
    ParametrizedSurfaceOutput&& a_rhs)
    : paraboloid_m(a_rhs.paraboloid_m),
      pt_from_bezier_split_m(std::move(a_rhs.pt_from_bezier_split_m)),
      arc_list_m(std::move(a_rhs.arc_list_m)) {}

inline ParametrizedSurfaceOutput& ParametrizedSurfaceOutput::operator=(
    ParametrizedSurfaceOutput&& a_rhs) {
  if (this != &a_rhs) {
    paraboloid_m = a_rhs.paraboloid_m;
    pt_from_bezier_split_m = std::move(a_rhs.pt_from_bezier_split_m);
    arc_list_m = std::move(a_rhs.arc_list_m);
  }
  return *this;
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

inline TriangulatedSurfaceOutput ParametrizedSurfaceOutput::triangulate(
    const double a_length_scale, const UnsignedIndex_t a_nsplit) const {
#ifndef IRL_USE_TRIANGLE
  std::cout << "IRL not compiled with the Triangle library. Exiting..."
            << std::endl;
  std::exit(-1);
#endif

  const UnsignedIndex_t nArcs = this->size();
  double length_scale = DBL_MAX;
  const auto& aligned_paraboloid = paraboloid_m.getAlignedParaboloid();

  std::cout << "Triangulating parametrized surface with " << nArcs
            << " arcs:" << std::endl;

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
      if (a_length_scale >= 0.0) {
        nSplit = static_cast<UnsignedIndex_t>(arc_length / a_length_scale);
        nSplit = nSplit < a_nsplit ? a_nsplit : nSplit;
      }
      const double step = 1.0 / static_cast<double>(nSplit);
      length_scale = std::min(length_scale, step * arc_length);
      if (a_length_scale >= 0.0) length_scale = a_length_scale;
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
    std::cout << "Passing " << input_points.size()
              << " points to the mesher, with maxarea = "
              << 0.5 * length_scale * length_scale << std::endl;
    in.numberofpoints = input_points.size() / 2;
    in.pointlist = input_points.data();

    std::vector<int> pointmarkerlist(in.numberofpoints, 1);
    in.pointmarkerlist = pointmarkerlist.data();

    in.numberofsegments = input_segments.size() / 2;
    in.segmentlist = input_segments.data();

    in.numberofholes = input_holes.size() / 2;
    if (in.numberofholes > 0) {
      in.holelist = input_holes.data();
    }

    char flags[50];
    sprintf(flags, "pzqYYia%.15fQ", 0.5 * length_scale * length_scale);
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
      const Pt base_pt = vertex + datum;
      vertex = Pt(0.0, 0.0, 0.0);
      for (UnsignedIndex_t d = 0; d < 3; ++d) {
        for (UnsignedIndex_t n = 0; n < 3; ++n) {
          vertex[n] += ref_frame[d][n] * base_pt[d];
        }
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
    free(in.segmentmarkerlist);
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
