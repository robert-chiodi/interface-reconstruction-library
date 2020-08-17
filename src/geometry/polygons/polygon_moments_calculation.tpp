// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_POLYGON_MOMENTS_CALCULATION_TPP_
#define SRC_GEOMETRY_POLYGONS_POLYGON_MOMENTS_CALCULATION_TPP_

#include <array>

namespace IRL {

template <class Derived, class VertexType, class SimplexType>
const Derived&
PolygonMomentsCalculationCommon<Derived, VertexType, SimplexType>::getDerived(
    void) const {
  return static_cast<const Derived&>(*this);
}

template <class Derived, class VertexType, class SimplexType>
UnsignedIndex_t PolygonMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::getNumberOfSimplicesInDecomposition(void)
    const {
  return this->getDerived().getNumberOfSimplicesInDecomposition();
}

template <class Derived, class VertexType, class SimplexType>
SimplexType PolygonMomentsCalculationCommon<Derived, VertexType, SimplexType>::
    getSimplexFromDecomposition(
        const UnsignedIndex_t a_tri_number_to_get) const {
  return this->getDerived().getSimplexFromDecomposition(a_tri_number_to_get);
}

template <class Derived, class VertexType, class SimplexType>
const VertexType&
    PolygonMomentsCalculationCommon<Derived, VertexType, SimplexType>::
    operator[](const UnsignedIndex_t a_index) const {
  return this->getDerived().access(a_index);
}

template <class Derived, class VertexType, class SimplexType>
UnsignedIndex_t PolygonMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::getNumberOfVertices(void) const {
  return this->getDerived().getNumberOfVerticesInObject();
}

template <class Derived, class VertexType, class SimplexType>
const Plane& PolygonMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::getPlaneOfExistence(void) const {
  return this->getDerived().getPlaneOfExistence_derived();
}

template <class Derived, class VertexType, class SimplexType>
inline Volume PolygonMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::calculateVolume(void) const {
  return IRL::calculateMoments((*this), Volume2D_Functor());
}

template <class Derived, class VertexType, class SimplexType>
inline VolumeMomentsAndNormal PolygonMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::calculateVolumeMomentsAndNormal(void)
    const {
  return IRL::calculateMoments((*this), VolumeMomentsAndNormal2D_Functor());
}

template <class Derived, class VertexType, class SimplexType>
inline Volume PolygonMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::calculateAbsoluteVolume(void) const {
  return std::fabs(this->calculateVolume());
}

template <class Derived, class VertexType, class SimplexType>
Volume PolygonMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::calculateConvexVolume(void) const {
  if (this->getNumberOfVertices() == 0) {
    return 0.0;
  }
  const Normal& plane_normal = this->getPlaneOfExistence().normal();
  auto dimensions_for_ascending_normal =
      this->getDimensionsOrderedForAscendingFaceNormalMagnitude(plane_normal);
  return std::fabs(this->calculate2DArea(dimensions_for_ascending_normal[0],
                                         dimensions_for_ascending_normal[1]) /
                   plane_normal[dimensions_for_ascending_normal[2]]);
}

template <class Derived, class VertexType, class SimplexType>
inline double PolygonMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::calculateSign(void) const {
  return std::copysign(1.0, this->calculateVolume());
}

template <class Derived, class VertexType, class SimplexType>
inline Pt PolygonMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::calculateCentroid(void) const {
  return IRL::calculateMoments((*this), Centroid2D_Functor());
}

template <class Derived, class VertexType, class SimplexType>
inline VolumeMoments PolygonMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::calculateMoments() const {
  return IRL::calculateMoments((*this), VolumeMoments2D_Functor());
}

template <class Derived, class VertexType, class SimplexType>
double PolygonMomentsCalculationCommon<Derived, VertexType, SimplexType>::
    calculate2DArea(const UnsignedIndex_t a_index_0,
                    const UnsignedIndex_t a_index_1) const {
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
std::array<UnsignedIndex_t, 3>
PolygonMomentsCalculationCommon<Derived, VertexType, SimplexType>::
    getDimensionsOrderedForAscendingFaceNormalMagnitude(Normal a_normal) const {
  std::array<UnsignedIndex_t, 3> indices{0, 1, 2};
  for (auto& element : a_normal) {
    element = std::fabs(element);
  }
  sortAscendingBasedOnOtherArray(&indices, &a_normal);
  return indices;
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYGONS_POLYGON_MOMENTS_CALCULATION_TPP_
