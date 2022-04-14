// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_TPP_
#define IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_TPP_

#include <float.h>
#include <cassert>
#include <cmath>

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/generic_cutting/paraboloid_intersection/wedge_computation.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/helpers/mymath.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

namespace IRL {

template <class ReturnType>
ReturnType computeV1Contribution(const Pt& a_ref_pt, const Pt& a_pt_0,
                                 const Pt& a_pt_1);

template <>
inline Volume computeV1Contribution(const Pt& a_ref_pt, const Pt& a_pt_0,
                                    const Pt& a_pt_1) {
  return (a_ref_pt[2] + a_pt_0[2] + a_pt_1[2]) / 6.0 *
         ((a_pt_0[0] - a_ref_pt[0]) * (a_pt_1[1] - a_ref_pt[1]) -
          (a_pt_1[0] - a_ref_pt[0]) * (a_pt_0[1] - a_ref_pt[1]));
}

template <>
inline VolumeMoments computeV1Contribution(const Pt& a_ref_pt, const Pt& a_pt_0,
                                           const Pt& a_pt_1) {
  std::cout << "Not yet implemented" << std::endl;
  std::exit(-1);
}

template <class ReturnType>
ReturnType computeV2Contribution(const AlignedParaboloid& a_aligned_paraboloid,
                                 const Pt& a_pt_0, const Pt& a_pt_1);

template <>
inline Volume computeV2Contribution(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1) {
  return (a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) / 12.0 *
         (-a_pt_0[2] - a_pt_1[2] +
          a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] +
          a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1]);
}

template <>
inline VolumeMoments computeV2Contribution(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_pt_0,
    const Pt& a_pt_1) {
  std::cout << "Not yet implemented" << std::endl;
  std::exit(-1);
}

// Starts from an entry, returns the exit that is reached.
template <class ReturnType, class HalfEdgeType>
ReturnType computeUnclippedSegmentV1Contribution(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_ref_pt,
    const HalfEdgeType a_entry_half_edge, HalfEdgeType& a_exit_half_edge,
    const bool skip_first = false) {
  ReturnType full_moments = ReturnType::fromScalarConstant(0.0);

  assert(a_entry_half_edge->getPreviousVertex()->isClipped() ||
         (a_entry_half_edge->getPreviousVertex()->doesNotNeedToSeek() &&
          a_entry_half_edge->getNextHalfEdge()->getVertex()->isNotClipped()));

  auto current_half_edge = a_entry_half_edge->getNextHalfEdge();
  auto vertex = current_half_edge->getVertex();
  if (skip_first && vertex->needsToSeek()) {
    current_half_edge = current_half_edge->getNextHalfEdge();
  }
  auto prev_pt = current_half_edge->getPreviousVertex()->getLocation().getPt();

  if (!(skip_first && vertex->doesNotNeedToSeek())) {
    while (true) {
      vertex = current_half_edge->getVertex();
      const auto& curr_pt = vertex->getLocation().getPt();
      full_moments +=
          computeV1Contribution<ReturnType>(a_ref_pt, prev_pt, curr_pt);
      if (vertex->needsToSeek()) {
        prev_pt = curr_pt;
        current_half_edge = current_half_edge->getNextHalfEdge();
      } else {
        assert(!(current_half_edge->getPreviousVertex()->isClipped() ||
                 (current_half_edge->getPreviousVertex()->doesNotNeedToSeek() &&
                  current_half_edge->getNextHalfEdge()
                      ->getVertex()
                      ->isNotClipped())));
        break;
      }
    }
  }
  a_exit_half_edge = current_half_edge;
  return full_moments;
}

template <class ReturnType, class HalfEdgeType, class SurfaceOutputType>
ReturnType computeNewEdgeSegmentContribution(
    const AlignedParaboloid& a_aligned_paraboloid, const Pt& a_ref_pt,
    const HalfEdgeType a_entry_half_edge, const HalfEdgeType a_exit_half_edge,
    const bool skip_first, SurfaceOutputType* a_surface) {
  ReturnType full_moments = ReturnType::fromScalarConstant(0.0);
  // Handle new edge on exit->entry
  if (!skip_first) {
    full_moments += computeV1Contribution<ReturnType>(
        a_ref_pt, a_exit_half_edge->getVertex()->getLocation().getPt(),
        a_entry_half_edge->getVertex()->getLocation().getPt());
  }

  full_moments += computeV2Contribution<ReturnType>(
      a_aligned_paraboloid,
      a_exit_half_edge->getVertex()->getLocation().getPt(),
      a_entry_half_edge->getVertex()->getLocation().getPt());

  full_moments += orientAndApplyWedgeCorrection<ReturnType>(
      a_aligned_paraboloid, a_exit_half_edge, a_entry_half_edge, a_surface);

  return full_moments;
}

template <class ReturnType, class HalfEdgeType, class SurfaceOutputType>
ReturnType computeContribution(const AlignedParaboloid& a_aligned_paraboloid,
                               const Pt& a_ref_pt,
                               const HalfEdgeType a_entry_half_edge,
                               const HalfEdgeType a_exit_half_edge,
                               const bool skip_first,
                               SurfaceOutputType* a_surface) {
  ReturnType full_moments = ReturnType::fromScalarConstant(0.0);

  full_moments += computeUnclippedSegmentV1Contribution(
      a_aligned_paraboloid, a_ref_pt, a_entry_half_edge, a_exit_half_edge,
      skip_first);
  full_moments += computeNewEdgeSegmentContribution(
      a_aligned_paraboloid, a_ref_pt, a_entry_half_edge, a_exit_half_edge,
      skip_first, a_surface);

  return full_moments;
}

template <class ReturnType>
ReturnType integrateFaceOnlyIntersection(const AlignedParaboloid& a_paraboloid,
                                         const Plane& a_face_plane);

template <>
inline Volume integrateFaceOnlyIntersection(
    const AlignedParaboloid& a_paraboloid, const Plane& a_face_plane) {
  assert(a_paraboloid.a() * a_paraboloid.b() > 0.0);
  assert(std::fabs(a_face_plane.normal()[2]) > DBL_EPSILON);
  const double a = -a_face_plane.normal()[0] / a_face_plane.normal()[2];
  const double b = -a_face_plane.normal()[1] / a_face_plane.normal()[2];
  const double c = a_face_plane.distance() / a_face_plane.normal()[2];
  const double factor = 4.0 * a_paraboloid.a() * a_paraboloid.b() * c -
                        a_paraboloid.a() * b * b - a_paraboloid.b() * a * a;
  return std::copysign(
      M_PI * factor * factor /
          (32.0 * std::pow(a_paraboloid.a() * a_paraboloid.b(), 2.5)),
      -a_face_plane.normal()[2]);
}

template <>
inline VolumeMoments integrateFaceOnlyIntersection(
    const AlignedParaboloid& a_paraboloid, const Plane& a_face_plane) {
  assert(a_paraboloid.a() * a_paraboloid.b() > 0.0);
  assert(std::fabs(a_face_plane.normal()[2]) > DBL_EPSILON);
  const double a = -a_face_plane.normal()[0] / a_face_plane.normal()[2];
  const double b = -a_face_plane.normal()[1] / a_face_plane.normal()[2];
  const double c = a_face_plane.distance() / a_face_plane.normal()[2];

  const double x_center = a / (2.0 * a_paraboloid.a());
  const double y_center = b / (2.0 * a_paraboloid.b());
  const double K = (4.0 * a_paraboloid.a() * a_paraboloid.b() * c +
                    a_paraboloid.a() * b * b + a_paraboloid.b() * a * a);

  auto moments = VolumeMoments::fromScalarConstant(0.0);
  moments.volume() = std::copysign(
      M_PI * fabs(K) * K /
          (32.0 * std::pow(a_paraboloid.a() * a_paraboloid.b(), 2.5)),
      -a_face_plane.normal()[2]);
  moments.centroid()[0] = x_center;
  moments.centroid()[1] = y_center;
  moments.centroid()[2] =
      (8.0 * a_paraboloid.a() * a_paraboloid.b() * c +
       5.0 * a_paraboloid.a() * b * b + 5.0 * a_paraboloid.b() * a * a) /
      (12.0 * a_paraboloid.a() * a_paraboloid.b());
  moments.multiplyByVolume();
  return moments;
}

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithParaboloid(SegmentedHalfEdgePolyhedronType* a_polytope,
                                  HalfEdgePolytopeType* a_complete_polytope,
                                  const Paraboloid& a_paraboloid) {
  // Move into reconstruction reference frame
  const UnsignedIndex_t original_number_of_vertices =
      a_polytope->getNumberOfVertices();
  const auto& datum = a_paraboloid.getDatum();
  const auto& ref_frame = a_paraboloid.getReferenceFrame();
  assert(ref_frame.isOrthonormalBasis());

  for (UnsignedIndex_t v = 0; v < original_number_of_vertices; ++v) {
    const auto& original_pt = a_polytope->getVertex(v)->getLocation().getPt();
    typename SegmentedHalfEdgePolyhedronType::pt_type projected_location;
    auto& pt = projected_location.getPt();
    for (UnsignedIndex_t n = 0; n < 3; ++n) {
      pt[n] = ref_frame[n] * original_pt - datum[n];
    }
    a_polytope->getVertex(v)->setLocation(projected_location);
  }

  // Recalculate in face plane information
  for (auto& face : (*a_polytope)) {
    auto normal = Normal(0.0, 0.0, 0.0);
    const auto starting_half_edge = face->getStartingHalfEdge();
    auto current_half_edge = starting_half_edge;
    auto next_half_edge = starting_half_edge->getNextHalfEdge();
    const auto& start_location =
        starting_half_edge->getPreviousVertex()->getLocation();
    do {
      normal += crossProduct(
          current_half_edge->getVertex()->getLocation() - start_location,
          next_half_edge->getVertex()->getLocation() - start_location);
      current_half_edge = next_half_edge;
      next_half_edge = next_half_edge->getNextHalfEdge();
    } while (next_half_edge != starting_half_edge);
    normal.normalize();
    face->setPlane(Plane(normal, normal * start_location));
  }

  // Compute intersection
  ReturnType moments;
  if constexpr (has_paraboloid_surface<ReturnType>::value) {
    moments.getSurface().setParaboloid(a_paraboloid);
    moments.getMoments() =
        intersectPolyhedronWithParaboloid<typename ReturnType::moment_type>(
            a_polytope, a_complete_polytope,
            a_paraboloid.getAlignedParaboloid(), &moments.getSurface());
  } else {
    NoSurfaceOutput* surf = nullptr;
    moments = intersectPolyhedronWithParaboloid<ReturnType>(
        a_polytope, a_complete_polytope, a_paraboloid.getAlignedParaboloid(),
        surf);
  }

  // Rotate base polyhedron back
  const UnsignedIndex_t number_of_vertices = a_polytope->getNumberOfVertices();
  for (UnsignedIndex_t v = 0; v < number_of_vertices; ++v) {
    const auto original_pt =
        a_polytope->getVertex(v)->getLocation().getPt() + datum;
    typename SegmentedHalfEdgePolyhedronType::pt_type projected_location;
    auto& pt = projected_location.getPt();
    pt = Pt(0.0, 0.0, 0.0);
    for (UnsignedIndex_t d = 0; d < 3; ++d) {
      for (UnsignedIndex_t n = 0; n < 3; ++n) {
        pt[n] += ref_frame[d][n] * original_pt[d];
      }
    }
    a_polytope->getVertex(v)->setLocation(projected_location);
  }

  return moments;
}

template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class SurfaceOutputType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
intersectPolyhedronWithParaboloid(SegmentedHalfEdgePolyhedronType* a_polytope,
                                  HalfEdgePolytopeType* a_complete_polytope,
                                  const AlignedParaboloid& a_paraboloid,
                                  SurfaceOutputType* a_surface) {
  // Assumed a_polytope is already rotated to be in same coordinate system
  // as a_paraboloid.

  // Below function computes the entire integration
  const auto moments = formParaboloidIntersectionBases<ReturnType>(
      a_polytope, a_complete_polytope, a_paraboloid, 0, a_surface);

  return moments;
}

template <class PtType, class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType>
void placeSingleIntercept(const PtType& a_intersection_location,
                          HalfEdgeType* a_half_edge_with_intersection,
                          SegmentedHalfEdgePolytopeType* a_polytope,
                          HalfEdgePolytopeType* a_complete_polytope) {
  auto first_intersection_vertex = a_complete_polytope->getNewVertex(
      typename SegmentedHalfEdgePolytopeType::vertex_type(
          a_intersection_location));
  first_intersection_vertex->setToSeek();
  a_polytope->addVertex(first_intersection_vertex);
  HalfEdgeType* new_half_edge = separateIntersectedHalfEdge(
      first_intersection_vertex, a_half_edge_with_intersection, a_polytope,
      a_complete_polytope);
  createOppositeHalfEdgeFromIntersection(a_half_edge_with_intersection,
                                         new_half_edge, a_complete_polytope);
}

template <class HalfEdgeType, class SegmentedHalfEdgePolytopeType,
          class HalfEdgePolytopeType, class VertexType>
void placeDoubleIntercept(const StackVector<std::pair<VertexType, double>, 2>&
                              a_intersection_location,
                          HalfEdgeType* a_half_edge_with_intersection,
                          SegmentedHalfEdgePolytopeType* a_polytope,
                          HalfEdgePolytopeType* a_complete_polytope) {
  assert(a_intersection_location.size() == 2);

  // Need to place furthest vertex first so that a_half_edge_with_intersection
  // remains attached to same vertex it started on. Also want to keep property
  // that the half edge the vertex stores has the previousVertex also
  // unclipped. This means for the first (furthest) vertex, need to reference
  // the opposite half edge of the current one. This enables more efficient
  // face truncation in the next phase of the algorithm, since new vertices
  // will always have a half edge going from unclipped->new->clipped.

  placeSingleIntercept(a_intersection_location[0].first,
                       a_half_edge_with_intersection, a_polytope,
                       a_complete_polytope);
  a_half_edge_with_intersection->getPreviousVertex()->setHalfEdge(
      a_half_edge_with_intersection->getOppositeHalfEdge());
  placeSingleIntercept(a_intersection_location[1].first,
                       a_half_edge_with_intersection, a_polytope,
                       a_complete_polytope);
}

template <class VertexType>
void checkAndFindIntercepts(
    const AlignedParaboloid& a_paraboloid, const VertexType& a_pt_0,
    const VertexType& a_pt_1,
    StackVector<std::pair<VertexType, double>, 2>* a_intercepts,
    const double a_nudge_epsilon = 0.0) {
  const auto pt_diff = a_pt_1 - a_pt_0;
  const double a = a_paraboloid.a() * pt_diff[0] * pt_diff[0] +
                   a_paraboloid.b() * pt_diff[1] * pt_diff[1];
  const double b = 2.0 * (a_paraboloid.a() * pt_diff[0] * a_pt_0[0] +
                          a_paraboloid.b() * pt_diff[1] * a_pt_0[1]) +
                   pt_diff[2];
  const double c = a_paraboloid.a() * a_pt_0[0] * a_pt_0[0] +
                   a_paraboloid.b() * a_pt_0[1] * a_pt_0[1] + a_pt_0[2];

  a_intercepts->resize(0);
  const auto solutions = solveQuadratic(a, b, c);
  for (const auto& solution : solutions) {
    if (solution >= -a_nudge_epsilon && solution <= 1.0 + a_nudge_epsilon) {
      a_intercepts->push_back(std::pair<VertexType, double>(
          VertexType(a_pt_0 + solution * pt_diff), solution));
    }
  }
}

template <class VertexType>
bool vertexBelow(const VertexType& a_pt,
                 const AlignedParaboloid& a_paraboloid) {
  return a_pt[2] < -(a_paraboloid.a() * a_pt[0] * a_pt[0] +
                     a_paraboloid.b() * a_pt[1] * a_pt[1]);
}

// This algorithm is based on the method shared from the website below, and
// from several stackoverflow posts citing that website.
// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
inline bool isPtBeforeIntersectionWithEdge(
    const std::array<double, 2>& a_test_pt, const Pt& a_vertex_0,
    const Pt& a_vertex_1) {
  if ((a_test_pt[1] > a_vertex_0[1]) == (a_test_pt[1] > a_vertex_1[1])) {
    return false;  // Projected ray never intersects edge.
  }
  double location_of_intersection_along_ray =
      (a_vertex_0[0] - a_vertex_1[0]) * (a_test_pt[1] - a_vertex_1[1]) /
          (a_vertex_0[1] - a_vertex_1[1]) +
      a_vertex_1[0];
  // Intersection was to the right if location_of_intersection_along_ray is
  // greater.
  return a_test_pt[0] < location_of_intersection_along_ray;
}

// If centroid is outside of polygon, or any vertex on face
// is inside the ellipse, then the ellipse is not contained by the face.
template <class HalfEdgeType>
bool ellipseContainedInFace(const AlignedParaboloid& a_aligned_paraboloid,
                            const Plane& a_face_plane,
                            HalfEdgeType* const a_half_edge) {
  const auto& face_normal = a_face_plane.normal();
  const std::array<double, 2> conic_center{
      {face_normal[0] / (2.0 * a_aligned_paraboloid.a() * face_normal[2]),
       face_normal[1] / (2.0 * a_aligned_paraboloid.b() * face_normal[2])}};
  const double delta_face = a_face_plane.distance() / face_normal[2];
  const double gamma_face =
      a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
      a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] - delta_face;
  if (a_aligned_paraboloid.a() * gamma_face < 0.0) {
    return false;
  }

  // First we will check if centroid is in the bounding box
  // of the face polygon. Due to the fact we are assuming
  // the paraboloid was in Z direction, we assume the ellipse
  // lives on the x/y plane and project the polygon
  // down to it as well (essentially neglecting the z component).
  auto current_half_edge = a_half_edge;
  std::array<double, 2> xy_min{{DBL_MAX, DBL_MAX}};
  std::array<double, 2> xy_max{{-DBL_MAX, -DBL_MAX}};
  do {
    const Pt& location = current_half_edge->getVertex()->getLocation().getPt();
    for (UnsignedIndex_t d = 0; d < 2; ++d) {
      xy_min[d] = std::min(xy_min[d], location[d]);
      xy_max[d] = std::max(xy_max[d], location[d]);
    }
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != a_half_edge);
  if (conic_center[0] < xy_min[0] || conic_center[0] > xy_max[0] ||
      conic_center[1] < xy_min[1] || conic_center[1] > xy_max[1]) {
    return false;
  }

  // If not outside bounding box, need to use
  // a more refined test based on the number
  // of edges a ray emitted from the point crosses.
  // Even edges means the point is outside the face.
  current_half_edge = a_half_edge;
  bool pt_internal_to_polygon = false;
  do {
    const Pt& location_0 =
        current_half_edge->getPreviousVertex()->getLocation().getPt();
    const Pt& location_1 =
        current_half_edge->getVertex()->getLocation().getPt();
    if (isPtBeforeIntersectionWithEdge(conic_center, location_0, location_1)) {
      pt_internal_to_polygon = !pt_internal_to_polygon;
    }
    current_half_edge = current_half_edge->getNextHalfEdge();
  } while (current_half_edge != a_half_edge);
  return pt_internal_to_polygon;
}

template <class ReturnType, class HalfEdgeType, class SurfaceOutputType>
ReturnType orientAndApplyWedgeCorrection(const AlignedParaboloid& a_paraboloid,
                                         HalfEdgeType* a_start,
                                         HalfEdgeType* a_end,
                                         SurfaceOutputType* a_surface) {
  ReturnType moments = ReturnType::fromScalarConstant(0.0);
  const auto& first_vertex_location =
      a_start->getVertex()->getLocation().getPt();
  const auto& second_vertex_location =
      a_end->getVertex()->getLocation().getPt();

  if (squaredMagnitude(second_vertex_location - first_vertex_location) <
      100.0 * DBL_EPSILON * DBL_EPSILON) {
    return moments;
  }

  auto face_plane = a_end->getFace()->getPlane();
  std::array<Normal, 2> tangents = {
      computeTangentVectorAtPoint(a_paraboloid, face_plane,
                                  first_vertex_location),
      computeTangentVectorAtPoint(a_paraboloid, face_plane,
                                  second_vertex_location)};

  const bool first_valid_edge_found =
      orientInitialTangents(a_start, face_plane, &tangents[0]);
  if (first_valid_edge_found) {
    const bool second_valid_edge_found =
        orientInitialTangents(a_end, face_plane, &tangents[1]);
    if (second_valid_edge_found) {
      moments = -(computeWedgeCorrection<ReturnType>(
          a_paraboloid, face_plane, first_vertex_location,
          second_vertex_location, tangents[0], tangents[1], a_surface));
    }
  }
  return moments;
}

template <class HalfEdgeType>
Normal determineNudgeDirection(const AlignedParaboloid& a_paraboloid,
                               HalfEdgeType* a_current_edge,
                               const Pt& a_inter_pt) {
  Normal nudge_direction = Normal(2.0 * a_paraboloid.a() * a_inter_pt[0],
                                  2.0 * a_paraboloid.b() * a_inter_pt[1], 1.0);
  Normal avg_face_normal =
      a_current_edge->getFace()->getPlane().normal() +
      a_current_edge->getOppositeHalfEdge()->getFace()->getPlane().normal();
  if (nudge_direction * avg_face_normal > 0.0) nudge_direction *= -1.0;
  nudge_direction.normalize();
  return nudge_direction;
}

// Note: This essentially abandons the half edges and vertices
// in a_complete_polytope, doing nothing to reclaim their memory usage.
// This should not be a problem unless resetPolyhedron is called
// many times for the same polyhedron intersection operation
template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, void>
resetPolyhedron(SegmentedHalfEdgePolyhedronType* a_polytope,
                HalfEdgePolytopeType* a_complete_polytope) {
  UnsignedIndex_t original_verts = 0;
  for (UnsignedIndex_t v = 0; v < a_polytope->getNumberOfVertices(); ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    if (vertex.needsToSeek()) {
      ++original_verts;
    } else {
      // New vertex from intersection, remove it and patch half-edges
      auto current_edge = vertex.getHalfEdge();
      auto original_edge = current_edge->getNextHalfEdge();
      doubleLinkHalfEdges(current_edge->getPreviousHalfEdge(), original_edge);
      original_edge->getFace()->setStartingHalfEdge(original_edge);
      auto opposite_edge = current_edge->getOppositeHalfEdge();
      doubleLinkHalfEdges(
          original_edge->getOppositeHalfEdge()->getPreviousHalfEdge(),
          opposite_edge);
      opposite_edge->getFace()->setStartingHalfEdge(opposite_edge);
      setMutualOpposites(original_edge, opposite_edge);
    }
  }
  a_polytope->setNumberOfVertices(original_verts);
  for (UnsignedIndex_t f = 0; f < a_polytope->getNumberOfFaces(); ++f) {
    (*a_polytope)[f]->clearIntersections();
  }
}

template <class SegmentedHalfEdgePolyhedronType, class HalfEdgePolytopeType>
void nudgePolyhedron(SegmentedHalfEdgePolyhedronType* a_polytope,
                     HalfEdgePolytopeType* a_complete_polytope,
                     const Normal a_nudge_direction,
                     const double a_nudge_epsilon,
                     const UnsignedIndex_t a_nudge_iter) {
  std::cout << " -----> Nudge along the direction " << a_nudge_direction
            << std::endl;
  for (UnsignedIndex_t n = 0; n < a_polytope->getNumberOfVertices(); ++n) {
    Pt shifted_location = a_polytope->getVertex(n)->getLocation() +
                          std::pow(1.1, static_cast<double>(a_nudge_iter)) *
                              a_nudge_epsilon * a_nudge_direction;
    a_polytope->getVertex(n)->setLocation(shifted_location);
  }
}

// Assumes paraboloid of function 0 = a*x^2 + b*y^2 + z.
// We will truncate the polyhedron to exist in the region
// q < a*x^2 + b*y^2 + z
template <class ReturnType, class SegmentedHalfEdgePolyhedronType,
          class HalfEdgePolytopeType, class SurfaceOutputType>
enable_if_t<is_polyhedron<SegmentedHalfEdgePolyhedronType>::value, ReturnType>
formParaboloidIntersectionBases(SegmentedHalfEdgePolyhedronType* a_polytope,
                                HalfEdgePolytopeType* a_complete_polytope,
                                const AlignedParaboloid& a_aligned_paraboloid,
                                const UnsignedIndex_t a_nudge_iter,
                                SurfaceOutputType* a_surface) {
  using vertex_type = typename SegmentedHalfEdgePolyhedronType::vertex_type;
  using face_type = typename SegmentedHalfEdgePolyhedronType::face_type;
  using half_edge_type = typename HalfEdgePolytopeType::half_edge_type;

  assert(!(a_surface != nullptr &&
           std::is_same<SurfaceOutputType, NoSurfaceOutput>::value));

  ReturnType full_moments = ReturnType::fromScalarConstant(0.0);

  // Initialising variables for handling degenerate cases
  bool requires_nudge = false;
  const double nudge_epsilon = 10.0 * DBL_EPSILON;
  Normal nudge_direction;

  // Identify elliptic case
  const bool elliptic =
      a_aligned_paraboloid.a() * a_aligned_paraboloid.b() > 0.0;

  // Mark vertices as clipped/unclipped
  const auto starting_number_of_vertices = a_polytope->getNumberOfVertices();
  const auto starting_number_of_faces = a_polytope->getNumberOfFaces();
  UnsignedIndex_t number_of_vertices_above = 0;
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    vertex.setAsUnnecessaryToSeek();  // Reset all
    if (vertexBelow(vertex.getLocation(), a_aligned_paraboloid)) {
      vertex.markToBeNotClipped();
    } else {
      vertex.markToBeClipped();
      ++number_of_vertices_above;
    }
  }

  // Early termination cases, only possible with elliptic
  if (elliptic && a_aligned_paraboloid.a() > 0.0 &&
      number_of_vertices_above == 0) {
    // Whole volume below
    return ReturnType::calculateMoments(a_polytope);
  }

  if (elliptic && a_aligned_paraboloid.a() < 0.0 &&
      number_of_vertices_above == starting_number_of_vertices) {
    // Zero volume - will be current value of full_moments
    return full_moments;
  }

  // Clear visitation knowledge from polytope.
  for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
    (*a_polytope)[f]->markAsNotVisited();
    (*a_polytope)[f]->clearIntersections();
  }

  // Will have 0, 1, or 2 intercepts.
  // If intercept exists, place into HalfEdgeStructure
  const bool check_from_clipped = true;
  //      !elliptic || (elliptic && a_aligned_paraboloid.a() > 0.0);
  const bool check_from_unclipped = true;
  //      !elliptic || (elliptic && !check_from_clipped);

  StackVector<std::pair<Pt, double>, 2> edge_intercepts;
  for (UnsignedIndex_t v = 0; v < starting_number_of_vertices; ++v) {
    auto& vertex = *(a_polytope->getVertex(v));
    vertex.setToSeek();
    if (check_from_clipped && vertex.isClipped()) {
      // CASE WHERE STARTING VERTEX IS CLIPPED
      auto current_edge = vertex.getHalfEdge();
      const auto starting_edge = current_edge;
      do {
        // If it has needsToSeek set, it means it is a new vertex or
        // already visited. Either way, do not need to check for intersection
        if (current_edge->getPreviousVertex()->needsToSeek()) {
          current_edge =
              current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
          continue;
        }
        // If previous vertex is not clipped, single-intercept
        // If previous vertex is clipped, need to check for double-intercept
        // Checking for double-intercept and calculating single-intercept
        // is the same routine, so just always do it.
        checkAndFindIntercepts(
            a_aligned_paraboloid,
            current_edge->getPreviousVertex()->getLocation().getPt(),
            current_edge->getVertex()->getLocation().getPt(), &edge_intercepts,
            nudge_epsilon);

        // assert(current_edge->getPreviousVertex()->isClipped() ||
        //        edge_intercepts.size() == 1);
        // Size of returned intercepts indicates single or double intercept
        // (or none)
        if (edge_intercepts.size() == 1) {
          // Check for intersection near end point
          if (edge_intercepts[0].second < nudge_epsilon ||
              (1.0 - edge_intercepts[0].second) < nudge_epsilon) {
            requires_nudge = true;
            nudge_direction = determineNudgeDirection(
                a_aligned_paraboloid, current_edge, edge_intercepts[0].first);
            break;
          }

          // Add vertex to half edge structure, resetting connectivity and
          // creating a new half edge (and new opposite half edge)
          placeSingleIntercept(edge_intercepts[0].first, current_edge,
                               a_polytope, a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->addIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->setStartingHalfEdge(opposite_half_edge);
          opposite_face->addIntersection();

        } else if (edge_intercepts.size() == 2) {
          // Check for intersection near end point
          if (edge_intercepts[0].second < nudge_epsilon) {
            requires_nudge = true;
            nudge_direction = determineNudgeDirection(
                a_aligned_paraboloid, current_edge, edge_intercepts[0].first);
            break;
          }
          if ((1.0 - edge_intercepts[1].second) < nudge_epsilon) {
            requires_nudge = true;
            nudge_direction = determineNudgeDirection(
                a_aligned_paraboloid, current_edge, edge_intercepts[1].first);
            break;
          }

          placeDoubleIntercept(edge_intercepts, current_edge, a_polytope,
                               a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->setStartingHalfEdge(
              current_edge->getPreviousHalfEdge()->getPreviousHalfEdge());
          current_face->addDoubleIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->setStartingHalfEdge(
              current_edge->getOppositeHalfEdge());
          opposite_face->addDoubleIntersection();
        }
        current_edge =
            current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
      } while (current_edge != starting_edge);
    } else if (check_from_unclipped && vertex.isNotClipped()) {
      // CASE WHERE STARTING VERTEX IS UNCLIPPED
      auto current_edge = vertex.getHalfEdge();
      const auto starting_edge = current_edge;
      do {
        // If it has needsToSeek set, it means it is a new vertex or
        // already visited. Either way, do not need to check for intersection
        if (current_edge->getPreviousVertex()->needsToSeek()) {
          current_edge =
              current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
          continue;
        }
        // If previous vertex is clipped, single-intercept
        // If previous vertex is not clipped, need to check for
        // double-intercept Checking for double-intercept and calculating
        // single-intercept is the same routine, so just always do it.
        checkAndFindIntercepts(
            a_aligned_paraboloid,
            current_edge->getPreviousVertex()->getLocation().getPt(),
            current_edge->getVertex()->getLocation().getPt(), &edge_intercepts,
            nudge_epsilon);

        // assert(current_edge->getPreviousVertex()->isNotClipped() ||
        //        edge_intercepts.size() == 1);
        // Size of returned intercepts indicates single or double intercept
        // (or none)
        if (edge_intercepts.size() == 1) {
          // Check for intersection near end point
          if (edge_intercepts[0].second < nudge_epsilon ||
              (1.0 - edge_intercepts[0].second) < nudge_epsilon) {
            requires_nudge = true;
            nudge_direction = determineNudgeDirection(
                a_aligned_paraboloid, current_edge, edge_intercepts[0].first);
            break;
          }

          // Add vertex to half edge structure, resetting connectivity and
          // creating a new half edge (and new opposite half edge)
          placeSingleIntercept(edge_intercepts[0].first, current_edge,
                               a_polytope, a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->setStartingHalfEdge(
              current_edge->getPreviousHalfEdge());
          current_face->addIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->addIntersection();

        } else if (edge_intercepts.size() == 2) {
          // Check for intersection near end point
          if (edge_intercepts[0].second < nudge_epsilon) {
            requires_nudge = true;
            nudge_direction = determineNudgeDirection(
                a_aligned_paraboloid, current_edge, edge_intercepts[0].first);
            break;
          }
          if ((1.0 - edge_intercepts[1].second) < nudge_epsilon) {
            requires_nudge = true;
            nudge_direction = determineNudgeDirection(
                a_aligned_paraboloid, current_edge, edge_intercepts[1].first);
            break;
          }

          placeDoubleIntercept(edge_intercepts, current_edge, a_polytope,
                               a_complete_polytope);

          auto current_face = current_edge->getFace();
          current_face->markAsVisited();
          current_face->setStartingHalfEdge(
              current_edge->getPreviousHalfEdge());
          current_face->addDoubleIntersection();

          auto opposite_half_edge = current_edge->getOppositeHalfEdge();
          auto opposite_face = opposite_half_edge->getFace();
          opposite_face->markAsVisited();
          opposite_face->setStartingHalfEdge(
              opposite_half_edge->getNextHalfEdge());
          opposite_face->addDoubleIntersection();
        }
        current_edge =
            current_edge->getOppositeHalfEdge()->getPreviousHalfEdge();
      } while (current_edge != starting_edge);
    }
    if (requires_nudge) {
      for (UnsignedIndex_t i = v + 1; i < starting_number_of_vertices; ++i) {
        a_polytope->getVertex(i)->setToSeek();
      }
      break;
    }
  }
  assert(a_polytope->checkValidHalfEdgeStructure());
  // After leaving above loop, all faces that were intersected have a
  // starting half edge that ends on a new (intersection) vertex and the
  // previous vertex is unclipped (or new). (i.e., the edge exists in the
  // unclipped portion.
  // All new vertices also have their half edge being one that has its
  // edge in the unclipped portion.

  const auto vertices_after_intersection = a_polytope->getNumberOfVertices();
  const auto new_intersection_vertices =
      vertices_after_intersection - starting_number_of_vertices;

  for (UnsignedIndex_t v = starting_number_of_vertices;
       v < vertices_after_intersection; ++v) {
    // Original vertices will be set as needsToSeek()
    // Reset newly created vertices will have doesNotNeedToSeek()
    a_polytope->getVertex(v)->setAsUnnecessaryToSeek();
  }

  // If intersection is too close to corners, shift polyhedron and try again
  if (requires_nudge) {
    // Clean half-edge structure by removing intersections
    if (new_intersection_vertices > 0) {
      resetPolyhedron(a_polytope, a_complete_polytope);
    }
    nudgePolyhedron(a_polytope, a_complete_polytope, nudge_direction,
                    nudge_epsilon, a_nudge_iter + 1);
    return formParaboloidIntersectionBases<ReturnType>(
        a_polytope, a_complete_polytope, a_aligned_paraboloid, a_nudge_iter + 1,
        a_surface);
  }

  // Check for face-only intersections. Can only happen with elliptic
  if (elliptic) {
    if (a_aligned_paraboloid.a() > 0.0) {
      for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
        auto& face = *(*a_polytope)[f];
        if (face.getNumberOfIntersections() > 0) {
          continue;
        }
        // Face will be FaceOnly intersect if
        // 1 - Any vertex on the face is below the paraboloid
        // 2 - No z-component on the face is below the maximum for
        // a downward opening elliptic paraboloid, which is 0.
        bool face_valid = false;
        const auto starting_half_edge = face.getStartingHalfEdge();
        auto current_half_edge = starting_half_edge;
        do {
          const auto& vertex = *(current_half_edge->getVertex());
          if (vertex.isNotClipped()) {
            face_valid = false;
            break;
          }
          if (vertex.getLocation().getPt()[2] < 0.0) {
            face_valid = true;
          }
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != starting_half_edge);
        if (!face_valid) {
          continue;
        }

        // Made it this far, check if there is a face-only intersection
        const auto& face_plane = face.getPlane();
        if (std::fabs(face_plane.normal()[2]) > DBL_EPSILON) {
          // Get ellipse on this face
          if (ellipseContainedInFace(a_aligned_paraboloid, face_plane,
                                     starting_half_edge)) {
            const auto face_integ = integrateFaceOnlyIntersection<ReturnType>(
                a_aligned_paraboloid, face_plane);
            full_moments += face_integ;
            // Return surface parametrization
            if constexpr (!std::is_same<SurfaceOutputType,
                                        NoSurfaceOutput>::value) {
              addEllipseToSurfaceOutput(a_aligned_paraboloid, face_plane,
                                        a_surface);
            }
          }
        }
      }
    } else {
      // Case where a_aligned_paraboloid.a() < 0.0)
      for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
        auto& face = *(*a_polytope)[f];
        if (face.getNumberOfIntersections() > 0) {
          continue;
        }
        // Face will be FaceOnly intersect if
        // 1 - Any vertex on the face is above the paraboloid
        // 2 - No z-component on the face is above the minimum for
        // an upward opening elliptic paraboloid, which is 0.
        bool face_valid = false;
        const auto starting_half_edge = face.getStartingHalfEdge();
        auto current_half_edge = starting_half_edge;
        do {
          const auto& vertex = *(current_half_edge->getVertex());
          if (vertex.isClipped()) {
            face_valid = false;
            break;
          }
          if (vertex.getLocation().getPt()[2] > 0.0) {
            face_valid = true;
          }
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != starting_half_edge);
        if (!face_valid) {
          continue;
        }

        // Made it this far, check if there is a face-only intersection
        const auto& face_plane = face.getPlane();
        if (std::fabs(face_plane.normal()[2]) > DBL_EPSILON) {
          // Get ellipse on this face
          if (ellipseContainedInFace(a_aligned_paraboloid, face_plane,
                                     starting_half_edge)) {
            full_moments += integrateFaceOnlyIntersection<ReturnType>(
                a_aligned_paraboloid, face_plane);
            if constexpr (!std::is_same<SurfaceOutputType,
                                        NoSurfaceOutput>::value) {
              addEllipseToSurfaceOutput(a_aligned_paraboloid, face_plane,
                                        a_surface);
            }
          }
        }
      }
    }
  }

  if (new_intersection_vertices == 0) {
    if (number_of_vertices_above == starting_number_of_vertices) {
      // All points above
      return full_moments;
    } else {
      // All points below
      return ReturnType::calculateMoments(a_polytope) + full_moments;
    }
  }

  for (UnsignedIndex_t f = 0; f < starting_number_of_faces; ++f) {
    auto& face = *(*a_polytope)[f];
    const auto& face_normal = face.getPlane().normal();
    auto starting_half_edge = face.getStartingHalfEdge();

    const auto intersection_size = face.getNumberOfIntersections();
    if (intersection_size == 0) {
      if (starting_half_edge->getVertex()->isNotClipped()) {
        // The face is entirely below
        const auto& ref_pt =
            starting_half_edge->getVertex()->getLocation().getPt();
        auto current_half_edge =
            starting_half_edge->getNextHalfEdge()->getNextHalfEdge();
        auto prev_pt =
            current_half_edge->getPreviousVertex()->getLocation().getPt();
        do {
          const auto& curr_pt =
              current_half_edge->getVertex()->getLocation().getPt();
          full_moments +=
              computeV1Contribution<ReturnType>(ref_pt, prev_pt, curr_pt);
          prev_pt = curr_pt;
          current_half_edge = current_half_edge->getNextHalfEdge();
        } while (current_half_edge != starting_half_edge);
      }
    } else if (intersection_size == 2) {
      const auto& ref_pt =
          starting_half_edge->getVertex()->getLocation().getPt();
      half_edge_type* exit_half_edge;
      full_moments += computeUnclippedSegmentV1Contribution<ReturnType>(
          a_aligned_paraboloid, ref_pt, starting_half_edge, exit_half_edge,
          true);
      full_moments += computeNewEdgeSegmentContribution<ReturnType>(
          a_aligned_paraboloid, ref_pt, starting_half_edge, exit_half_edge,
          true, a_surface);
    } else {
      // More than 2 intersections. More complicated cases.
      const bool elliptic_face =
          a_aligned_paraboloid.a() * a_aligned_paraboloid.b() > 0.0 &&
          face_normal[2] != 0.0;
      const bool hyperbolic_face =
          a_aligned_paraboloid.a() * a_aligned_paraboloid.b() < 0.0 &&
          face_normal[2] != 0.0;

      // -----> If NO surface ouput : DOES NOT NEED SORTING
      if constexpr (std::is_same<SurfaceOutputType, NoSurfaceOutput>::value) {
        const bool reverse = a_aligned_paraboloid.a() < 0.0;
        if (elliptic_face) {
          const auto& ref_pt =
              starting_half_edge->getVertex()->getLocation().getPt();
          half_edge_type* exit_half_edge;
          auto current_edge = starting_half_edge;
          bool skip_first = true;
          std::size_t found_intersections = 0;
          do {
            full_moments += computeUnclippedSegmentV1Contribution<ReturnType>(
                a_aligned_paraboloid, ref_pt, current_edge, exit_half_edge,
                skip_first);

            if (reverse) {
              full_moments += computeNewEdgeSegmentContribution<ReturnType>(
                  a_aligned_paraboloid, ref_pt, current_edge, exit_half_edge,
                  skip_first, a_surface);
              current_edge = exit_half_edge->getNextHalfEdge();
              while (current_edge->getVertex()->needsToSeek()) {
                current_edge = current_edge->getNextHalfEdge();
              }
            } else {
              current_edge = exit_half_edge->getNextHalfEdge();
              while (current_edge->getVertex()->needsToSeek()) {
                current_edge = current_edge->getNextHalfEdge();
              }
              full_moments += computeNewEdgeSegmentContribution<ReturnType>(
                  a_aligned_paraboloid, ref_pt, current_edge, exit_half_edge,
                  current_edge == starting_half_edge, a_surface);
            }
            skip_first = false;
            found_intersections += 2;
          } while (found_intersections != intersection_size);
        } else {
          SmallVector<half_edge_type*, 10> intersections;
          // Find intersections and determine status
          auto current_edge = starting_half_edge;
          bool reverse_order = false;
          half_edge_type* exit_half_edge;
          std::size_t found_intersections = 0;
          const auto& ref_pt =
              starting_half_edge->getVertex()->getLocation().getPt();
          bool skip_first = true;
          do {
            full_moments += computeUnclippedSegmentV1Contribution<ReturnType>(
                a_aligned_paraboloid, ref_pt, current_edge, exit_half_edge,
                skip_first);
            current_edge->getVertex()->markAsEntry();
            exit_half_edge->getVertex()->markAsExit();
            intersections.push_back(current_edge);
            intersections.push_back(exit_half_edge);
            current_edge = exit_half_edge->getNextHalfEdge();
            while (current_edge->getVertex()->needsToSeek()) {
              current_edge = current_edge->getNextHalfEdge();
            }
            found_intersections += 2;
            skip_first = false;
          } while (found_intersections != intersection_size);

          if (hyperbolic_face) {
            const std::array<double, 2> conic_center{
                {face_normal[0] /
                     (2.0 * a_aligned_paraboloid.a() * face_normal[2]),
                 face_normal[1] /
                     (2.0 * a_aligned_paraboloid.b() * face_normal[2])}};
            const auto& pt_in_plane =
                intersections[0]->getVertex()->getLocation().getPt();
            const double delta_face = (face_normal[0] * pt_in_plane[0] +
                                       face_normal[1] * pt_in_plane[1] +
                                       face_normal[2] * pt_in_plane[2]) /
                                      face_normal[2];
            const double gamma_face =
                a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
                a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] -
                delta_face;
            const double z_diff =
                face_normal[0] * conic_center[0] / face_normal[2] +
                face_normal[1] * conic_center[1] / face_normal[2] - gamma_face -
                2.0 * delta_face;
            const std::size_t split_ind =
                a_aligned_paraboloid.a() * gamma_face > 0.0 ? 0 : 1;

            bool vertices_on_same_branch = true;
            for (std::size_t i = 1; i < intersection_size; ++i) {
              const auto& curr_pt =
                  intersections[i]->getVertex()->getLocation().getPt();
              if ((pt_in_plane[split_ind] - conic_center[split_ind]) *
                      (curr_pt[split_ind] - conic_center[split_ind]) <
                  0.0) {
                vertices_on_same_branch = false;
                break;
              }
            }
            bool total_invert = z_diff < 0.0 ? !vertices_on_same_branch
                                             : vertices_on_same_branch;
            if (total_invert) {
              reverse_order = true;
            }
          } else {
            Pt intersection_avg = Pt(0.0, 0.0, 0.0);
            for (std::size_t i = 0; i < intersection_size; ++i) {
              const auto& curr_pt =
                  intersections[i]->getVertex()->getLocation().getPt();
              intersection_avg += curr_pt;
            }
            intersection_avg /= static_cast<double>(intersection_size);
            const double z_diff = intersection_avg[2] +
                                  a_aligned_paraboloid.a() *
                                      intersection_avg[0] *
                                      intersection_avg[0] +
                                  a_aligned_paraboloid.b() *
                                      intersection_avg[1] * intersection_avg[1];
            if (z_diff > 0.0) {
              reverse_order = true;
            }
          }

          // Traverse face from entry->exit
          // Integrate internal portion of face
          // Integrate wedge on new edge from exit->entry
          bool first_entry = true;
          // Loop over vertices, only use entry vertices
          for (std::size_t i = 0; i < intersection_size; i += 2) {
            // Identify entry vertex
            auto entry_half_edge = intersections[i];
            // Loop over internal portion from entry->exit
            int exit_index;
            if (reverse_order) {
              exit_index = i != intersection_size - 1 ? i + 1 : 0;
            } else {
              exit_index = i != 0 ? i - 1 : intersection_size - 1;
            }

            auto exit_half_edge = intersections[exit_index];
            full_moments += computeNewEdgeSegmentContribution<ReturnType>(
                a_aligned_paraboloid, ref_pt, entry_half_edge, exit_half_edge,
                first_entry, a_surface);
            first_entry = false;
          }
          // Clear list of intersections
          intersections.clear();
        }
      }
      // -----> If surface ouput : NEEDS SORTING
      else {
        using stype = std::pair<half_edge_type*, double>;
        SmallVector<stype, 10> intersections;
        // Find intersections and determine status
        auto current_edge = starting_half_edge;
        bool reverse_order = false;
        half_edge_type* exit_half_edge;
        std::size_t found_intersections = 0;
        const auto& ref_pt =
            starting_half_edge->getVertex()->getLocation().getPt();
        bool skip_first = true;
        do {
          full_moments += computeUnclippedSegmentV1Contribution<ReturnType>(
              a_aligned_paraboloid, ref_pt, current_edge, exit_half_edge,
              skip_first);
          current_edge->getVertex()->markAsEntry();
          exit_half_edge->getVertex()->markAsExit();
          intersections.push_back(
              std::pair<half_edge_type*, double>({current_edge, DBL_MAX}));
          intersections.push_back(
              std::pair<half_edge_type*, double>({exit_half_edge, DBL_MAX}));
          current_edge = exit_half_edge->getNextHalfEdge();
          while (current_edge->getVertex()->needsToSeek()) {
            current_edge = current_edge->getNextHalfEdge();
          }
          found_intersections += 2;
          skip_first = false;
        } while (found_intersections != intersection_size);

        if (elliptic_face) {
          const std::array<double, 2> conic_center{
              {face_normal[0] /
                   (2.0 * a_aligned_paraboloid.a() * face_normal[2]),
               face_normal[1] /
                   (2.0 * a_aligned_paraboloid.b() * face_normal[2])}};
          const double normal_invert = std::copysign(1.0, face_normal[2]);
          const double invert =
              a_aligned_paraboloid.a() < 0.0 ? -normal_invert : normal_invert;

          for (auto& element : intersections) {
            const auto& pt = element.first->getVertex()->getLocation().getPt();
            element.second = invert * std::atan2(pt[1] - conic_center[1],
                                                 pt[0] - conic_center[0]);
          }
          std::sort(intersections.begin(), intersections.end(),
                    [](const stype& a, const stype& b) {
                      return a.second < b.second;
                    });
        } else if (hyperbolic_face) {
          std::size_t pos_end = 0;
          std::size_t neg_end = 0;
          auto intersection_copy = intersections;
          const std::array<double, 2> conic_center{
              {face_normal[0] /
                   (2.0 * a_aligned_paraboloid.a() * face_normal[2]),
               face_normal[1] /
                   (2.0 * a_aligned_paraboloid.b() * face_normal[2])}};
          const auto& pt_in_plane =
              intersections[0].first->getVertex()->getLocation().getPt();
          const double delta_face = (face_normal[0] * pt_in_plane[0] +
                                     face_normal[1] * pt_in_plane[1] +
                                     face_normal[2] * pt_in_plane[2]) /
                                    face_normal[2];
          const double gamma_face =
              a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] +
              a_aligned_paraboloid.b() * conic_center[1] * conic_center[1] -
              delta_face;
          const std::size_t split_ind =
              a_aligned_paraboloid.a() * gamma_face > 0.0 ? 0 : 1;
          const std::size_t store_ind = split_ind == 0 ? 1 : 0;
          const double z_center_plane =
              -face_normal[0] * conic_center[0] / face_normal[2] -
              face_normal[1] * conic_center[1] / face_normal[2] + delta_face;
          const double z_center_paraboloid =
              -a_aligned_paraboloid.a() * conic_center[0] * conic_center[0] -
              a_aligned_paraboloid.b() * conic_center[1] * conic_center[1];
          double total_invert = std::copysign(1.0, face_normal[2]);
          total_invert = z_center_plane < z_center_paraboloid ? total_invert
                                                              : -total_invert;
          total_invert = split_ind == 0 ? total_invert : -total_invert;
          for (auto& element : intersections) {
            const auto& pt = element.first->getVertex()->getLocation().getPt();
            if (pt[split_ind] > conic_center[split_ind]) {
              const auto loc = pos_end++;
              intersection_copy[loc].first = element.first;
              intersection_copy[loc].second = total_invert * pt[store_ind];
            } else {
              const auto loc = intersection_size - 1 - (neg_end++);
              intersection_copy[loc].first = element.first;
              intersection_copy[loc].second = -total_invert * pt[store_ind];
            }
          }
          assert(pos_end + neg_end == intersection_size);
          std::sort(intersection_copy.begin(),
                    intersection_copy.begin() + pos_end,
                    [](const stype& a, const stype& b) {
                      return a.second < b.second;
                    });
          std::sort(intersection_copy.begin() + pos_end,
                    intersection_copy.end(),
                    [](const stype& a, const stype& b) {
                      return a.second < b.second;
                    });
          intersections = intersection_copy;
        } else {
          Pt intersection_avg = Pt(0.0, 0.0, 0.0);
          for (std::size_t i = 0; i < intersection_size; ++i) {
            const auto& curr_pt =
                intersections[i].first->getVertex()->getLocation().getPt();
            intersection_avg += curr_pt;
          }
          intersection_avg /= static_cast<double>(intersection_size);
          const double z_diff = intersection_avg[2] +
                                a_aligned_paraboloid.a() * intersection_avg[0] *
                                    intersection_avg[0] +
                                a_aligned_paraboloid.b() * intersection_avg[1] *
                                    intersection_avg[1];
          if (face_normal[2] != 0.0) {
            const double normal_invert = std::copysign(1.0, face_normal[2]);
            const double invert = z_diff < 0.0 ? normal_invert : -normal_invert;
            for (auto& element : intersections) {
              const auto& pt =
                  element.first->getVertex()->getLocation().getPt();
              element.second = invert * std::atan2(pt[1] - intersection_avg[1],
                                                   pt[0] - intersection_avg[0]);
            }
            std::sort(intersections.begin(), intersections.end(),
                      [](const stype& a, const stype& b) {
                        return a.second < b.second;
                      });
          } else {
          }
        }
        // Traverse face from entry->exit
        // Integrate internal portion of face
        // Integrate wedge on new edge from exit->entry
        // Loop over vertices, only use entry vertices
        auto prev_vertex = intersections[0].first->getPreviousVertex();
        auto entry_half_edge = intersections[0].first;
        const bool entry_first =
            (prev_vertex->isClipped() ||
             (prev_vertex->doesNotNeedToSeek() &&
              entry_half_edge->getNextHalfEdge()->getVertex()->isNotClipped()));
        std::size_t start_id = entry_first ? 0 : 1;
        for (std::size_t i = start_id; i < intersection_size; i += 2) {
          // Identify entry vertex
          const auto entry_half_edge = intersections[i].first;
          // Loop over internal portion from entry->exit
          const UnsignedIndex_t exit_index =
              i != 0 ? i - 1 : intersection_size - 1;
          const auto exit_half_edge = intersections[exit_index].first;
          full_moments += computeNewEdgeSegmentContribution<ReturnType>(
              a_aligned_paraboloid, ref_pt, entry_half_edge, exit_half_edge,
              false, a_surface);
        }
        // Clear list of intersections
        intersections.clear();
      }
    }
  }

  return full_moments;
}
}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_PARABOLOID_INTERSECTION_TPP_
