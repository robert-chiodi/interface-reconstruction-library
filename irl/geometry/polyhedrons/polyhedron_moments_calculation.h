// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_POLYHEDRONS_POLYHEDRON_MOMENTS_CALCULATION_H_
#define IRL_GEOMETRY_POLYHEDRONS_POLYHEDRON_MOMENTS_CALCULATION_H_

#include "irl/geometry/general/moment_calculation_through_simplices.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/pt_with_data.h"
#include "irl/moments/general_moments.h"
#include "irl/moments/volume.h"
#include "irl/moments/volume_moments.h"
#include "irl/moments/volume_moments_and_doubles.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class GeometryType>
class ProxyTet;

template <class Derived, class VertexType, class SimplexType>
class PolyhedronMomentsCalculationCommon {
  const Derived& getDerived(void) const;

 public:
  UnsignedIndex_t getNumberOfSimplicesInDecomposition(void) const;

  SimplexType getSimplexFromDecomposition(
      const UnsignedIndex_t a_tet_number_to_get) const;

  /// \brief Calculate and return volume of the tet.
  inline Volume calculateVolume(void) const;

  /// \brief Calculate and return signed volume of the tet.
  /// See Owkes & Desjardins, JCP, 2014.
  inline Volume calculateAbsoluteVolume(void) const;

  /// \brief Calculate sign for the tet
  /// See Owkes & Desjardins, JCP, 2014.
  inline double calculateSign(void) const;

  /// \brief Calculate and return centroid of the tet.
  inline Pt calculateCentroid(void) const;

  /// \brief Calculate and return volume weighted VolumeMoments.
  inline VolumeMoments calculateMoments() const;

  /// \brief Calculate General Moments
  template <std::size_t ORDER>
  inline GeneralMoments3D<ORDER> calculateGeneralMoments(void) const;
};

template <class Derived, class VertexType, class SimplexType>
class PolyhedronMomentsCalculation
    : public PolyhedronMomentsCalculationCommon<Derived, VertexType,
                                                SimplexType> {};

template <class Derived, class FunctorType, UnsignedIndex_t kArrayLength,
          class SimplexType>
class PolyhedronMomentsCalculation<
    Derived, PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>,
    SimplexType>
    : public PolyhedronMomentsCalculationCommon<
          Derived, PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>,
          SimplexType> {
 public:
  inline VolumeMomentsAndDoubles<kArrayLength> calculateVolumeMomentsAndDoubles(
      void) const;
};

}  // namespace IRL

#include "irl/geometry/polyhedrons/polyhedron_moments_calculation.tpp"

#endif  // IRL_GEOMETRY_POLYHEDRONS_POLYHEDRON_MOMENTS_CALCULATION_H_
