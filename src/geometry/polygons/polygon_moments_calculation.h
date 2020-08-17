// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYGONS_POLYGON_MOMENTS_CALCULATION_H_
#define SRC_GEOMETRY_POLYGONS_POLYGON_MOMENTS_CALCULATION_H_

#include "src/geometry/general/moment_calculation_through_simplices.h"
#include "src/geometry/general/pt.h"
#include "src/moments/volume.h"
#include "src/moments/volume_moments.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class Derived, class VertexType, class SimplexType>
class PolygonMomentsCalculationCommon {
  const Derived& getDerived(void) const;

 public:
  UnsignedIndex_t getNumberOfSimplicesInDecomposition(void) const;

  SimplexType getSimplexFromDecomposition(
      const UnsignedIndex_t a_tri_number_to_get) const;

  const VertexType& operator[](const UnsignedIndex_t a_index) const;

  UnsignedIndex_t getNumberOfVertices(void) const;

  const Plane& getPlaneOfExistence(void) const;

  /// \brief Calculate and return volume of the tri.
  inline Volume calculateVolume(void) const;

  inline Volume calculateAbsoluteVolume(void) const;

  inline Volume calculateConvexVolume(void) const;

  inline double calculateSign(void) const;

  /// \brief Calculate and return centroid of the tri.
  inline Pt calculateCentroid(void) const;

  /// \brief Calculate and return volume weighted VolumeMoments.
  inline VolumeMoments calculateMoments() const;

  /// \brief Calculate and return volume weighted VolumeMoments.
  inline VolumeMomentsAndNormal calculateVolumeMomentsAndNormal() const;

 private:
  double calculate2DArea(const UnsignedIndex_t a_index_0,
                         const UnsignedIndex_t a_index_1) const;

  std::array<UnsignedIndex_t, 3>
  getDimensionsOrderedForAscendingFaceNormalMagnitude(Normal a_normal) const;
};

template <class Derived, class VertexType, class SimplexType>
class PolygonMomentsCalculation
    : public PolygonMomentsCalculationCommon<Derived, VertexType, SimplexType> {
};

}  // namespace IRL

#include "src/geometry/polygons/polygon_moments_calculation.tpp"

#endif  // SRC_GEOMETRY_POLYGONS_POLYGON_MOMENTS_CALCULATION_H_
