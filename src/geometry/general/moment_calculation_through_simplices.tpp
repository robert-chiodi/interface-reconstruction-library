// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GEOMETRY_GENERAL_MOMENT_CALCULATION_THROUGH_SIMPLICES_TPP_
#define SRC_GEOMETRY_GENERAL_MOMENT_CALCULATION_THROUGH_SIMPLICES_TPP_

#include "src/geometry/general/normal.h"
#include "src/helpers/mymath.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class GeometryType, class CalculationFunctor>
auto calculateMoments(const GeometryType& a_geometry,
                      CalculationFunctor a_moment_accumulator) ->
    typename CalculationFunctor::ReturnType {
  for (UnsignedIndex_t s = 0;
       s < a_geometry.getNumberOfSimplicesInDecomposition(); ++s) {
    a_moment_accumulator(a_geometry.getSimplexFromDecomposition(s));
  }
  return a_moment_accumulator.getMoments();
}

template <class SimplexType>
void Volume3D_Functor::operator()(const SimplexType& a_simplex) {
  const auto& datum = a_simplex[3];
  volume_m += scalarTripleProduct(a_simplex[0].getPt() - datum.getPt(),
		                          a_simplex[1].getPt() - datum.getPt(),
				                  a_simplex[2].getPt() - datum.getPt());
}

Volume3D_Functor::ReturnType Volume3D_Functor::getMoments(void) const {
  return volume_m / 6.0;
}

template <class SimplexType>
void VolumeMoments3D_Functor::operator()(const SimplexType& a_simplex) {
  const auto& datum = a_simplex[3];
  auto six_times_tet_volume =
		  scalarTripleProduct(a_simplex[0].getPt() - datum.getPt(),
		  		              a_simplex[1].getPt() - datum.getPt(),
							  a_simplex[2].getPt() - datum.getPt());
  volume_moments_m.volume() += six_times_tet_volume;
  volume_moments_m.centroid() +=
      six_times_tet_volume * (a_simplex[0].getPt() + a_simplex[1].getPt() +
                              a_simplex[2].getPt() + datum.getPt());
}

VolumeMoments3D_Functor::ReturnType VolumeMoments3D_Functor::getMoments(
    void) const {
  auto moments_copy = volume_moments_m;
  moments_copy /= 6.0;
  moments_copy.centroid() *= 0.25;
  return moments_copy;
}

template <class SimplexType>
void Centroid3D_Functor::operator()(const SimplexType& a_simplex) {
  volume_moments_functor_m(a_simplex);
}

Centroid3D_Functor::ReturnType Centroid3D_Functor::getMoments(void) const {
  auto volume_moments = volume_moments_functor_m.getMoments();
  volume_moments.normalizeByVolume();
  return volume_moments.centroid();
}

template <UnsignedIndex_t kArrayLength>
template <class SimplexType>
void VolumeMomentsAndDoubles3D_Functor<kArrayLength>::operator()(
    const SimplexType& a_simplex) {
  const auto& datum = a_simplex[3];
  auto six_times_tet_volume =
       scalarTripleProduct(a_simplex[0].getPt() - datum.getPt(),
		                   a_simplex[1].getPt() - datum.getPt(),
		                   a_simplex[2].getPt() - datum.getPt());
  volume_moments_and_doubles_m.volume() += six_times_tet_volume;
  volume_moments_and_doubles_m.centroid() +=
      six_times_tet_volume * (a_simplex[0].getPt() + a_simplex[1].getPt() +
                              a_simplex[2].getPt() + datum.getPt());
  for (UnsignedIndex_t n = 0; n < volume_moments_and_doubles_m.data().size();
       ++n) {
    volume_moments_and_doubles_m.data()[n] +=
        six_times_tet_volume *
        (a_simplex[0].getData()[n] + a_simplex[1].getData()[n] +
         a_simplex[2].getData()[n] + datum.getData()[n]);
  }
}

template <UnsignedIndex_t kArrayLength>
typename VolumeMomentsAndDoubles3D_Functor<kArrayLength>::ReturnType
VolumeMomentsAndDoubles3D_Functor<kArrayLength>::getMoments(void) const {
  auto moments_copy = volume_moments_and_doubles_m;
  moments_copy /= 6.0;
  moments_copy.centroid() *= 0.25;
  for (UnsignedIndex_t n = 0; n < moments_copy.data().size(); ++n) {
    moments_copy.data()[n] *= 0.25;
  }
  return moments_copy;
}

template <class SimplexType>
void Volume2D_Functor::operator()(const SimplexType& a_simplex) {
  const auto& datum_pt = a_simplex[0];
  Normal edge_cross_product =
      Normal::fromPt(crossProduct(a_simplex[1].getPt() - datum_pt.getPt(),
                                  a_simplex[2].getPt() - datum_pt.getPt()));

  // Dot product imparts signed-ness to the polygon area.
  double twice_times_tri_volume = magnitude(edge_cross_product);
  edge_cross_product.normalize();
  volume_m +=
      twice_times_tri_volume *
      dotProduct(a_simplex.getPlaneOfExistence().normal(), edge_cross_product);
}

Volume2D_Functor::ReturnType Volume2D_Functor::getMoments(void) const {
  return 0.5 * volume_m;
}

template <class SimplexType>
void VolumeMoments2D_Functor::operator()(const SimplexType& a_simplex) {
  const auto& datum_pt = a_simplex[0];
  Normal edge_cross_product =
      Normal::fromPt(crossProduct(a_simplex[1].getPt() - datum_pt.getPt(),
                                  a_simplex[2].getPt() - datum_pt.getPt()));

  // Dot product imparts signed-ness to the polygon area.
  double twice_times_tri_volume = magnitude(edge_cross_product);
  edge_cross_product.normalize();
  twice_times_tri_volume *=
      dotProduct(a_simplex.getPlaneOfExistence().normal(), edge_cross_product);

  volume_moments_m.volume() += twice_times_tri_volume;
  volume_moments_m.centroid() +=
      twice_times_tri_volume *
      (datum_pt.getPt() + a_simplex[1].getPt() + a_simplex[2].getPt());
}

VolumeMoments2D_Functor::ReturnType VolumeMoments2D_Functor::getMoments(
    void) const {
  auto moments_copy = volume_moments_m;
  moments_copy *= 0.5;
  moments_copy.centroid() /= 3.0;
  return moments_copy;
}

template <class SimplexType>
void Centroid2D_Functor::operator()(const SimplexType& a_simplex) {
  volume_moments_functor_m(a_simplex);
}

Centroid2D_Functor::ReturnType Centroid2D_Functor::getMoments(void) const {
  auto volume_moments = volume_moments_functor_m.getMoments();
  volume_moments.normalizeByVolume();
  return volume_moments.centroid();
}

template <class SimplexType>
void VolumeMomentsAndNormal2D_Functor::operator()(
    const SimplexType& a_simplex) {
  volume_moments_functor_m(a_simplex);
  normal_m = a_simplex.getPlaneOfExistence().normal();
}

VolumeMomentsAndNormal2D_Functor::ReturnType
VolumeMomentsAndNormal2D_Functor::getMoments(void) const {
  auto volume_moments = volume_moments_functor_m.getMoments();
  return {volume_moments, volume_moments.volume() * normal_m};
}

}  // namespace IRL

#endif  // SRC_GEOMETRY_GENERAL_MOMENT_CALCULATION_THROUGH_SIMPLICES_TPP_
