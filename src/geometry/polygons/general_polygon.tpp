// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_GENERAL_POLYGON_TPP_
#define SRC_GEOMETRY_POLYGONS_GENERAL_POLYGON_TPP_

#include <float.h>
#include <cmath>
#include <numeric>

#include "src/geometry/polygons/tri.h"
#include "src/helpers/geometric_cutting_helpers.h"
#include "src/helpers/helper.h"
#include "src/helpers/mymath.h"

namespace IRL {

namespace polygon_detail {
template <class PtType>
Plane getPlaneFromEdge(const PtType& a_starting_point,
                       const PtType& a_ending_point,
                       const Normal& a_face_normal,
                       const Pt& a_surface_centroid);

template <class PtType>
Plane getPlaneFromEdge(const PtType& a_starting_point,
                       const PtType& a_ending_point,
                       const Normal& a_face_normal,
                       const Pt& a_surface_centroid) {
  const auto edge_normal =
      Normal::fromPtNormalized(a_ending_point - a_starting_point);
  auto plane_normal = crossProductNormalized(edge_normal, a_face_normal);
  makeNormalFaceOutwardsFromVolume(a_surface_centroid,
                                   0.5 * (a_starting_point + a_ending_point),
                                   &plane_normal);
  return Plane(plane_normal, plane_normal * a_starting_point);
}

inline bool validNormal(const Normal& a_normal) {
  double normal_magnitude = magnitude(a_normal);
  return normal_magnitude > 0.9 || normal_magnitude < DBL_MIN;
}
}  // namespace polygon_detail

template <class Derived, class VertexType, class SimplexType>
Derived& GeneralPolygon<Derived, VertexType, SimplexType>::getDerived(void) {
  return static_cast<Derived&>(*this);
}
template <class Derived, class VertexType, class SimplexType>
const Derived& GeneralPolygon<Derived, VertexType, SimplexType>::getDerived(
    void) const {
  return static_cast<const Derived&>(*this);
}

template <class Derived, class VertexType, class SimplexType>
VertexType& GeneralPolygon<Derived, VertexType, SimplexType>::operator[](
    const UnsignedIndex_t a_index) {
  return this->getDerived().access(a_index);
}

template <class Derived, class VertexType, class SimplexType>
const VertexType& GeneralPolygon<Derived, VertexType, SimplexType>::operator[](
    const UnsignedIndex_t a_index) const {
  return this->getDerived().access(a_index);
}

template <class Derived, class VertexType, class SimplexType>
UnsignedIndex_t
GeneralPolygon<Derived, VertexType,
               SimplexType>::getNumberOfSimplicesInDecomposition(void) const {
  return this->getDerived().getNumberOfSimplicesInDecomposition();
}

template <class Derived, class VertexType, class SimplexType>
std::array<UnsignedIndex_t, 3>
GeneralPolygon<Derived, VertexType, SimplexType>::
    getSimplexIndicesFromDecomposition(
        const UnsignedIndex_t a_tet_number_to_get) const {
  return this->getDerived().getSimplexIndicesFromDecomposition(
      a_tet_number_to_get);
}

template <class Derived, class VertexType, class SimplexType>
SimplexType
GeneralPolygon<Derived, VertexType, SimplexType>::getSimplexFromDecomposition(
    const UnsignedIndex_t a_tet_number_to_get) const {
  return this->getDerived().getSimplexFromDecomposition(a_tet_number_to_get);
}

template <class Derived, class VertexType, class SimplexType>
UnsignedIndex_t GeneralPolygon<Derived, VertexType,
                               SimplexType>::getNumberOfVertices(void) const {
  return this->getDerived().getNumberOfVerticesInObject();
}

template <class Derived, class VertexType, class SimplexType>
void GeneralPolygon<Derived, VertexType, SimplexType>::setNumberOfVertices(
    const UnsignedIndex_t a_size) {
  this->getDerived().setNumberOfVerticesInObject(a_size);
}

template <class Derived, class VertexType, class SimplexType>
IRL::Pt GeneralPolygon<Derived, VertexType, SimplexType>::getLowerLimits(
    void) const {
  Pt pt_to_return(DBL_MAX, DBL_MAX, DBL_MAX);
  for (const auto& pt : (*this)) {
    pt_to_return[0] = std::min(pt_to_return[0], pt[0]);
    pt_to_return[1] = std::min(pt_to_return[1], pt[1]);
    pt_to_return[2] = std::min(pt_to_return[2], pt[2]);
  }
  return pt_to_return;
}

template <class Derived, class VertexType, class SimplexType>
IRL::Pt GeneralPolygon<Derived, VertexType, SimplexType>::getUpperLimits(
    void) const {
  Pt pt_to_return(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for (const auto& pt : (*this)) {
    pt_to_return[0] = std::max(pt_to_return[0], pt[0]);
    pt_to_return[1] = std::max(pt_to_return[1], pt[1]);
    pt_to_return[2] = std::max(pt_to_return[2], pt[2]);
  }
  return pt_to_return;
}

template <class Derived, class VertexType, class SimplexType>
typename GeneralPolygon<Derived, VertexType, SimplexType>::iterator
GeneralPolygon<Derived, VertexType, SimplexType>::begin(void) noexcept {
  return iterator(this->getDerived(), 0);
}
template <class Derived, class VertexType, class SimplexType>
typename GeneralPolygon<Derived, VertexType, SimplexType>::const_iterator
GeneralPolygon<Derived, VertexType, SimplexType>::begin(void) const noexcept {
  return this->cbegin();
}
template <class Derived, class VertexType, class SimplexType>
typename GeneralPolygon<Derived, VertexType, SimplexType>::const_iterator
GeneralPolygon<Derived, VertexType, SimplexType>::cbegin(void) const noexcept {
  return const_iterator(this->getDerived(), 0);
}
template <class Derived, class VertexType, class SimplexType>
typename GeneralPolygon<Derived, VertexType, SimplexType>::iterator
GeneralPolygon<Derived, VertexType, SimplexType>::end(void) noexcept {
  return iterator(this->getDerived(), this->getNumberOfVertices());
}
template <class Derived, class VertexType, class SimplexType>
typename GeneralPolygon<Derived, VertexType, SimplexType>::const_iterator
GeneralPolygon<Derived, VertexType, SimplexType>::end(void) const noexcept {
  return this->cend();
}
template <class Derived, class VertexType, class SimplexType>
typename GeneralPolygon<Derived, VertexType, SimplexType>::const_iterator
GeneralPolygon<Derived, VertexType, SimplexType>::cend(void) const noexcept {
  return const_iterator(this->getDerived(), this->getNumberOfVertices());
}

template <class Derived, class VertexType, class SimplexType>
void GeneralPolygon<Derived, VertexType,
                    SimplexType>::calculateAndSetPlaneOfExistence(void) {
  Normal normal = this->calculateNormal();
  this->setPlaneOfExistence(Plane(normal, normal * (*this)[0]));
}

template <class Derived, class VertexType, class SimplexType>
void GeneralPolygon<Derived, VertexType, SimplexType>::setPlaneOfExistence(
    const Plane& a_plane) {
  this->getDerived().setPlaneOfExistence_derived(a_plane);
}

template <class Derived, class VertexType, class SimplexType>
const Plane& GeneralPolygon<Derived, VertexType,
                            SimplexType>::getPlaneOfExistence(void) const {
  return this->getDerived().getPlaneOfExistence_derived();
}

template <class Derived, class VertexType, class SimplexType>
double GeneralPolygon<Derived, VertexType, SimplexType>::calculateSign(
    void) const {
  assert(magnitude(this->calculateNormal()) < 1.0e-10 ||
         magnitude(this->getPlaneOfExistence().normal()) > 0.9999);
  assert(std::fabs(this->getPlaneOfExistence().normal() *
                   this->calculateNormal()) > 0.9999);
  return this->getPlaneOfExistence().normal() * this->calculateNormal();
}

template <class Derived, class VertexType, class SimplexType>
Normal GeneralPolygon<Derived, VertexType, SimplexType>::calculateNormal(
    void) const {
  assert(this->getNumberOfVertices() == 0 ||
         polygon_detail::validNormal(Normal::fromPtNormalized(
             crossProduct((*this)[1] - (*this)[0], (*this)[2] - (*this)[0]))));
  return this->getNumberOfVertices() > 0
             ? Normal::fromPtNormalized(crossProduct((*this)[1] - (*this)[0],
                                                     (*this)[2] - (*this)[0]))
             : Normal(0.0, 0.0, 0.0);
}

template <class Derived, class VertexType, class SimplexType>
std::array<UnsignedIndex_t, 3>
GeneralPolygon<Derived, VertexType, SimplexType>::
    getDimensionsOrderedForAscendingFaceNormalMagnitude(Normal a_normal) const {
  std::array<UnsignedIndex_t, 3> indices;
  std::iota(indices.begin(), indices.end(), 0);
  for (auto& element : a_normal) {
    element = std::fabs(element);
  }
  sortAscendingBasedOnOtherArray(&indices, &a_normal);
  return indices;
}

template <class Derived, class VertexType, class SimplexType>
double GeneralPolygon<Derived, VertexType, SimplexType>::calculate2DArea(
    const UnsignedIndex_t a_index_0, const UnsignedIndex_t a_index_1) const {
  assert(this->getNumberOfVertices() > 0);
  UnsignedIndex_t nvert = this->getNumberOfVertices();
  double area = {0.0};
  for (UnsignedIndex_t n = 0; n < nvert - 1; ++n) {
    area += (*this)[n][a_index_0] * (*this)[(n + 1)][a_index_1] -
            (*this)[n][a_index_1] * (*this)[(n + 1)][a_index_0];
  }
  area += (*this)[nvert - 1][a_index_0] * (*this)[0][a_index_1] -
          (*this)[nvert - 1][a_index_1] * (*this)[0][a_index_0];
  return 0.5 * area;
}

template <class Derived, class VertexType, class SimplexType>
Pt GeneralPolygon<Derived, VertexType, SimplexType>::
    calculateNearestPtOnSurface(const Pt& a_pt_to_project) const {
  assert(this->getNumberOfVertices() > 0);
  const auto polygon_normal = this->calculateNormal();
  Pt projected_point = this->getPtOnPlane(a_pt_to_project, polygon_normal);
  if (!isPtInsidePolygon(projected_point, polygon_normal)) {
    projected_point = this->getNearestProjectedPtOnEdges(a_pt_to_project);
  }
  return projected_point;
}

template <class Derived, class VertexType, class SimplexType>
PlanarLocalizer GeneralPolygon<Derived, VertexType, SimplexType>::getLocalizer(
    void) const {
  PlanarLocalizer localizer_to_return;
  localizer_to_return.setNumberOfPlanes(this->getNumberOfVertices());
  auto face_normal = this->calculateNormal();
  auto polygon_centroid = this->calculateCentroid();
  for (UnsignedIndex_t edge = 0; edge < this->getNumberOfVertices() - 1;
       ++edge) {
    localizer_to_return[edge] = polygon_detail::getPlaneFromEdge(
        (*this)[edge], (*this)[edge + 1], face_normal, polygon_centroid);
  }
  localizer_to_return[this->getNumberOfVertices() - 1] =
      polygon_detail::getPlaneFromEdge((*this)[this->getNumberOfVertices() - 1],
                                       (*this)[0], face_normal,
                                       polygon_centroid);
  return localizer_to_return;
}

template <class Derived, class VertexType, class SimplexType>
void GeneralPolygon<Derived, VertexType, SimplexType>::reversePtOrdering(void) {
  Derived old_polygon = static_cast<const Derived&>(*this);
  for (UnsignedIndex_t n = 1; n < this->getNumberOfVertices(); ++n) {
    (*this)[n] = old_polygon[this->getNumberOfVertices() - n];
  }
}

template <class Derived, class VertexType, class SimplexType>
Pt GeneralPolygon<Derived, VertexType, SimplexType>::getPtOnPlane(
    const Pt& a_pt_to_project, const Normal& a_polygon_normal) const {
  double normal_distance_to_plane =
      Pt(a_pt_to_project - (*this)[0]) * a_polygon_normal;
  return a_pt_to_project -
         (normal_distance_to_plane * Normal::toPt(a_polygon_normal));
}

template <class Derived, class VertexType, class SimplexType>
bool GeneralPolygon<Derived, VertexType, SimplexType>::isPtInsidePolygon(
    const Pt& a_pt, const Normal& a_polygon_normal) const {
  auto dimensions_to_use =
      this->getDimensionsOrderedForAscendingFaceNormalMagnitude(
          a_polygon_normal);
  return isPtInternalTo2DPolygon(a_pt, dimensions_to_use);
}

// This algorithm is based on the method shared from the website below, and from
// several stackoverflow posts citing that website.
// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
template <class Derived, class VertexType, class SimplexType>
bool GeneralPolygon<Derived, VertexType, SimplexType>::isPtInternalTo2DPolygon(
    const Pt& a_pt,
    const std::array<UnsignedIndex_t, 3> a_dimension_masking) const {
  bool pt_internal_to_polygon = false;
  // If Pt is in the bounding box, check via edge intersections
  if (this->isPtInBoundingBox(a_pt)) {
    for (UnsignedIndex_t edge = 0; edge < this->getNumberOfVertices() - 1;
         ++edge) {
      if (isPtBeforeIntersectionWithEdge(a_pt, (*this)[edge], (*this)[edge + 1],
                                         a_dimension_masking)) {
        // Works because if pt is internal, should have an odd number
        // of intersections with the polygon edges
        pt_internal_to_polygon = !pt_internal_to_polygon;
      }
    }
    if (isPtBeforeIntersectionWithEdge(a_pt,
                                       (*this)[this->getNumberOfVertices() - 1],
                                       (*this)[0], a_dimension_masking)) {
      pt_internal_to_polygon = !pt_internal_to_polygon;
    }
  }
  return pt_internal_to_polygon;
}

template <class Derived, class VertexType, class SimplexType>
bool GeneralPolygon<Derived, VertexType, SimplexType>::isPtInBoundingBox(
    const Pt& a_pt) const {
  Pt lower_bounding_pt = this->getLowerLimits();
  Pt upper_bounding_pt = this->getUpperLimits();
  return (a_pt[0] >= lower_bounding_pt[0] && a_pt[0] <= upper_bounding_pt[0] &&
          a_pt[1] >= lower_bounding_pt[1] && a_pt[1] <= upper_bounding_pt[1] &&
          a_pt[2] >= lower_bounding_pt[2] && a_pt[2] <= upper_bounding_pt[2]);
}

template <class Derived, class VertexType, class SimplexType>
bool GeneralPolygon<Derived, VertexType, SimplexType>::
    isPtBeforeIntersectionWithEdge(
        const Pt& a_test_pt, const Pt& a_vertex_0, const Pt& a_vertex_1,
        const std::array<UnsignedIndex_t, 3> a_dimension_masking) const {
  if ((a_test_pt[a_dimension_masking[1]] >
       a_vertex_0[a_dimension_masking[1]]) ==
      (a_test_pt[a_dimension_masking[1]] >
       a_vertex_1[a_dimension_masking[1]])) {
    return false;  // Projected ray never intersects edge.
  }
  double location_of_intersection_along_ray =
      (a_vertex_0[a_dimension_masking[0]] -
       a_vertex_1[a_dimension_masking[0]]) *
          (a_test_pt[a_dimension_masking[1]] -
           a_vertex_1[a_dimension_masking[1]]) /
          (a_vertex_0[a_dimension_masking[1]] -
           a_vertex_1[a_dimension_masking[1]]) +
      a_vertex_0[a_dimension_masking[0]];
  // Intersection was to the right if location_of_intersection_along_ray is
  // greater.
  return a_test_pt[a_dimension_masking[0]] < location_of_intersection_along_ray;
}

template <class Derived, class VertexType, class SimplexType>
Pt GeneralPolygon<Derived, VertexType, SimplexType>::
    getNearestProjectedPtOnEdges(const Pt& a_pt_to_project) const {
  Pt current_nearest_point;
  double current_shortest_squared_distance = DBL_MAX;
  for (UnsignedIndex_t edge = 0; edge < this->getNumberOfVertices() - 1;
       ++edge) {
    Pt point_on_edge = this->getProjectedPtOnEdge(
        a_pt_to_project, (*this)[edge], (*this)[edge + 1]);
    this->takePtIfCloserToOriginalPt(a_pt_to_project, &point_on_edge,
                                     &current_nearest_point,
                                     &current_shortest_squared_distance);
  }
  Pt point_on_edge = this->getProjectedPtOnEdge(
      a_pt_to_project, (*this)[this->getNumberOfVertices() - 1], (*this)[0]);
  this->takePtIfCloserToOriginalPt(a_pt_to_project, &point_on_edge,
                                   &current_nearest_point,
                                   &current_shortest_squared_distance);
  return current_nearest_point;
}

template <class Derived, class VertexType, class SimplexType>
Pt GeneralPolygon<Derived, VertexType, SimplexType>::getProjectedPtOnEdge(
    const Pt& a_pt_to_project, const Pt& a_edge_pt_0,
    const Pt& a_edge_pt_1) const {
  Vec3 edge_vector = (a_edge_pt_1 - a_edge_pt_0);
  double edge_length_squared = squaredMagnitude(edge_vector);
  if (edge_length_squared < DBL_MIN) {
    return a_edge_pt_0;
  }
  Vec3 vertex_to_pt_vector = (a_pt_to_project - a_edge_pt_0);
  double distance_along_line =
      (edge_vector * vertex_to_pt_vector) / edge_length_squared;
  distance_along_line = clipBetween(0.0, distance_along_line, 1.0);
  return a_edge_pt_0 + distance_along_line * Pt::fromVec3(edge_vector);
}

template <class Derived, class VertexType, class SimplexType>
void GeneralPolygon<Derived, VertexType, SimplexType>::
    takePtIfCloserToOriginalPt(
        const Pt& a_pt_to_project, Pt* a_point_from_the_polygon_edge,
        Pt* a_current_closest_pt,
        double* a_current_shortest_squared_distance) const {
  double squared_distance_from_projected_pt_to_original =
      squaredMagnitude(a_pt_to_project - *a_point_from_the_polygon_edge);
  if (squared_distance_from_projected_pt_to_original <
      *a_current_shortest_squared_distance) {
    *a_current_shortest_squared_distance =
        squared_distance_from_projected_pt_to_original;
    *a_current_closest_pt = *a_point_from_the_polygon_edge;
  }
}

template <class Derived, class VertexType, class SimplexType>
inline std::ostream& operator<<(
    std::ostream& out,
    const GeneralPolygon<Derived, VertexType, SimplexType>& a_polygon_base) {
  out << "Polygon has " << a_polygon_base.getNumberOfVertices()
      << "vertices \n";
  for (UnsignedIndex_t vertex = 0;
       vertex < a_polygon_base.getNumberOfVertices(); ++vertex) {
    out << "Vertex " << vertex << " : " << a_polygon_base[vertex] << '\n';
  }
  return out;
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYGONS_GENERAL_POLYGON_TPP_
