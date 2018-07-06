// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_POLYHEDRONS_POLYHEDRON_MOMENTS_CALCULATION_TPP_
#define SRC_GEOMETRY_POLYHEDRONS_POLYHEDRON_MOMENTS_CALCULATION_TPP_

namespace IRL {

template <class Derived, class VertexType, class SimplexType>
const Derived& PolyhedronMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::getDerived(void) const {
  return static_cast<const Derived&>(*this);
}

template <class Derived, class VertexType, class SimplexType>
UnsignedIndex_t PolyhedronMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::getNumberOfSimplicesInDecomposition(void)
    const {
  return this->getDerived().getNumberOfSimplicesInDecomposition();
}

template <class Derived, class VertexType, class SimplexType>
SimplexType
PolyhedronMomentsCalculationCommon<Derived, VertexType, SimplexType>::
    getSimplexFromDecomposition(
        const UnsignedIndex_t a_tet_number_to_get) const {
  return this->getDerived().getSimplexFromDecomposition(a_tet_number_to_get);
}

template <class Derived, class VertexType, class SimplexType>
inline Volume PolyhedronMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::calculateVolume(void) const {
  return IRL::calculateMoments((*this), Volume3D_Functor());
}

template <class Derived, class VertexType, class SimplexType>
inline Volume PolyhedronMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::calculateAbsoluteVolume(void) const {
  return std::fabs(this->calculateVolume());
}

template <class Derived, class VertexType, class SimplexType>
inline double PolyhedronMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::calculateSign(void) const {
  return std::copysign(1.0, this->calculateVolume());
}

template <class Derived, class VertexType, class SimplexType>
inline Pt PolyhedronMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::calculateCentroid(void) const {
  return IRL::calculateMoments((*this), Centroid3D_Functor());
}

template <class Derived, class VertexType, class SimplexType>
inline VolumeMoments PolyhedronMomentsCalculationCommon<
    Derived, VertexType, SimplexType>::calculateMoments(void) const {
  return IRL::calculateMoments((*this), VolumeMoments3D_Functor());
}

template <class Derived, class FunctorType, UnsignedIndex_t kArrayLength,
          class SimplexType>
inline VolumeMomentsAndDoubles<kArrayLength> PolyhedronMomentsCalculation<
    Derived, PtWithDoublesStatelessFunctor<FunctorType, kArrayLength>,
    SimplexType>::calculateVolumeMomentsAndDoubles(void) const {
  return IRL::calculateMoments(
      (*this), VolumeMomentsAndDoubles3D_Functor<kArrayLength>());
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_POLYHEDRONS_POLYHEDRON_MOMENTS_CALCULATION_TPP_
