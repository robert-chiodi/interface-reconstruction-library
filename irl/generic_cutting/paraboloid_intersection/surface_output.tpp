// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_SURFACE_OUTPUT_TPP_
#define IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_SURFACE_OUTPUT_TPP_

namespace IRL {
template <class ScalarType>
inline void addEllipseToSurfaceOutput(
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid,
    const PlaneBase<ScalarType>& a_face_plane,
    ParametrizedSurfaceOutput* a_surface) {
  const auto aligned_paraboloid =
      AlignedParaboloidBase<double>(a_aligned_paraboloid);
  const auto face_plane = PlaneBase<double>(
      NormalBase<double>(static_cast<double>(a_face_plane.normal()[0]),
                         static_cast<double>(a_face_plane.normal()[1]),
                         static_cast<double>(a_face_plane.normal()[2])),
      static_cast<double>(a_face_plane.distance()));

  assert(aligned_paraboloid.a() * aligned_paraboloid.b() > 0.0);
  assert(std::fabs(face_plane.normal()[2]) > DBL_EPSILON);
  const double alpha = aligned_paraboloid.a();
  const double beta = aligned_paraboloid.b();
  const double a = -face_plane.normal()[0] / face_plane.normal()[2];
  const double b = -face_plane.normal()[1] / face_plane.normal()[2];
  const double c = face_plane.distance() / face_plane.normal()[2];
  const double semi_axis_x =
      std::sqrt((a * a / (4.0 * alpha) + b * b / (4.0 * beta) - c) / alpha);
  const double semi_axis_y =
      std::sqrt((a * a / (4.0 * alpha) + b * b / (4.0 * beta) - c) / beta);
  const double enclosing_radius =
      std::sqrt(semi_axis_x * semi_axis_x + semi_axis_y * semi_axis_y);
  const double x_center = -a / (2.0 * aligned_paraboloid.a());
  const double y_center = -b / (2.0 * aligned_paraboloid.b());
  static constexpr UnsignedIndex_t nSplit = 3;  // Minimum acceptable is 3
  std::array<Pt, nSplit> points;
  const double angle_interval = 2.0 * M_PI / static_cast<double>(nSplit);
  const double normal_invert = std::copysign(1.0, face_plane.normal()[2]);
  const double invert =
      aligned_paraboloid.a() < 0.0 ? normal_invert : -normal_invert;
  for (UnsignedIndex_t i = 0; i < nSplit; ++i) {
    const double theta = invert * static_cast<double>(i) * angle_interval;
    points[i] = Pt(x_center + semi_axis_x * cos(theta),
                   y_center + semi_axis_y * sin(theta), 0.0);
    points[i][2] = a * points[i][0] + b * points[i][1] + c;
  }
  auto gradF_0 = getParaboloidSurfaceNormal(aligned_paraboloid, points[0]);
  Normal tangent_0 = crossProduct(face_plane.normal(), gradF_0);
  tangent_0.normalize();
  Normal start_tangent = tangent_0;
  for (UnsignedIndex_t i = 0; i < nSplit - 1; ++i) {
    auto gradF_end =
        getParaboloidSurfaceNormal(aligned_paraboloid, points[i + 1]);
    Normal end_tangent = crossProduct(face_plane.normal(), gradF_end);
    end_tangent.normalize();
    if (dotProduct(start_tangent, points[i + 1] - points[i]) < 0.0) {
      start_tangent *= -1.0;
    }
    if (dotProduct(end_tangent, points[i + 1] - points[i]) > 0.0) {
      end_tangent *= -1.0;
    }
    a_surface->addArc(RationalBezierArc(points[i], start_tangent, points[i + 1],
                                        end_tangent, face_plane.normal(),
                                        aligned_paraboloid));
    start_tangent = -end_tangent;
  }
  if (dotProduct(tangent_0, points[0] - points[nSplit - 1]) > 0.0) {
    tangent_0 *= -1.0;
  }
  a_surface->addArc(RationalBezierArc(points[nSplit - 1], start_tangent,
                                      points[0], tangent_0, face_plane.normal(),
                                      aligned_paraboloid));
}
}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_SURFACE_OUTPUT_TPP_
