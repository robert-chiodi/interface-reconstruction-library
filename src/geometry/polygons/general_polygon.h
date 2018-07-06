// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_GENERAL_POLYGON_H_
#define SRC_GEOMETRY_POLYGONS_GENERAL_POLYGON_H_

#include <algorithm>

#include "src/geometry/polygons/polygon_moments_calculation.h"
#include "src/planar_reconstruction/planar_localizer.h"

namespace IRL {

template <class ContainerType>
class IteratorThroughBracketOperator;

template <class ContainerType>
class ConstIteratorThroughBracketOperator;

template <class Derived, class VertexType, class SimplexType>
class GeneralPolygon
    : public PolygonMomentsCalculation<Derived, VertexType, SimplexType> {
  Derived& getDerived(void);
  const Derived& getDerived(void) const;

  using iterator = IteratorThroughBracketOperator<Derived>;
  using const_iterator = ConstIteratorThroughBracketOperator<Derived>;

 public:
  using pt_type = VertexType;
  using value_t = pt_type;

  VertexType& operator[](const UnsignedIndex_t a_index);

  const VertexType& operator[](const UnsignedIndex_t a_index) const;

  UnsignedIndex_t getNumberOfSimplicesInDecomposition(void) const;

  std::array<UnsignedIndex_t, 3> getSimplexIndicesFromDecomposition(
      const UnsignedIndex_t a_tet) const;

  SimplexType getSimplexFromDecomposition(
      const UnsignedIndex_t a_tet_number_to_get) const;

  UnsignedIndex_t getNumberOfVertices(void) const;

  void setNumberOfVertices(const UnsignedIndex_t a_size);

  /// \brief Return a point for the lower limits of the polyhedron in 3D space.
  IRL::Pt getLowerLimits(void) const;
  /// \brief Return a point for the upper limits of the polyhedron in 3D space.
  IRL::Pt getUpperLimits(void) const;
  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  void calculateAndSetPlaneOfExistence(void);

  void setPlaneOfExistence(const Plane& a_plane);
  const Plane& getPlaneOfExistence(void) const;

  double calculateSign(void) const;
  /// \brief Const version of calculateNormal.
  Normal calculateNormal(void) const;

  /// \brief Reverse point ordering.
  void reversePtOrdering(void);

  /// \brief Returns a planar reconstruction that is equivalent to the
  /// polygon extruded infinitely along it's face normal.
  PlanarLocalizer getLocalizer(void) const;

  /// \brief Return closest Pt  to `a_pt_to_project` that lays on the surface.
  Pt calculateNearestPtOnSurface(const Pt& a_pt_to_project) const;

 private:
  /// \brief Calculate 2D area of polygon projected into 2D plane
  /// with greatest magnitude normal component for the plane.
  double calculate2DArea(const UnsignedIndex_t a_index_0,
                         const UnsignedIndex_t a_index_1) const;

  std::array<UnsignedIndex_t, 3>
  getDimensionsOrderedForAscendingFaceNormalMagnitude(Normal a_normal) const;

  Pt getPtOnPlane(const Pt& a_pt_to_project,
                  const Normal& a_polygon_normal) const;

  bool isPtInsidePolygon(const Pt& a_pt, const Normal& a_polygon_normal) const;

  bool isPtInternalTo2DPolygon(
      const Pt& a_pt,
      const std::array<UnsignedIndex_t, 3> a_dimension_masking) const;

  bool isPtInBoundingBox(const Pt& a_pt) const;

  bool isPtBeforeIntersectionWithEdge(
      const Pt& a_test_pt, const Pt& a_vertex_0, const Pt& a_vertex_1,
      const std::array<UnsignedIndex_t, 3> a_dimension_masking) const;

  Pt getNearestProjectedPtOnEdges(const Pt& a_pt_to_project) const;

  Pt getProjectedPtOnEdge(const Pt& a_pt_to_project, const Pt& a_edge_pt_0,
                          const Pt& a_edge_pt_1) const;

  void removePt(const UnsignedIndex_t a_index) {
    this->getDerived().removePt(a_index);
  }

  void takePtIfCloserToOriginalPt(
      const Pt& a_pt_to_project, Pt* a_point_from_the_polygon_edge,
      Pt* a_current_closest_pt,
      double* a_current_shortest_squared_distance) const;
};

template <class Derived, class VertexType, class SimplexType>
inline std::ostream& operator<<(
    std::ostream& out,
    const GeneralPolygon<Derived, VertexType, SimplexType>& a_polygon_base);

}  // namespace IRL

#include "src/geometry/polygons/general_polygon.tpp"

#endif  // SRC_GEOMETRY_POLYGONS_GENERAL_POLYGON_H_
