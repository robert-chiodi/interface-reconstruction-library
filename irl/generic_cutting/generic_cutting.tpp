// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_GENERIC_CUTTING_TPP_
#define IRL_GENERIC_CUTTING_GENERIC_CUTTING_TPP_

#include "irl/generic_cutting/analytic/rectangular_cuboid.h"
#include "irl/generic_cutting/analytic/tet.h"
#include "irl/moments/volume.h"
#include "irl/planar_reconstruction/planar_separator.h"

namespace IRL {

namespace generic_cutting_details {

template <class EncompassingType>
static constexpr enable_if_t<is_polyhedron<EncompassingType>::value, bool>
polytopeIsValid(const EncompassingType& a_polytope) {
  return true;
}

template <class EncompassingType>
enable_if_t<is_polygon<EncompassingType>::value, bool> polytopeIsValid(
    const EncompassingType& a_polytope) {
  return a_polytope.calculateAbsoluteVolume() <
             global_constants::MINIMUM_SURFACE_AREA_TO_TRACK ||
         a_polytope.getPlaneOfExistence().normal().calculateMagnitude() > 0.999;
}

template <class EncompassingType>
static constexpr enable_if_t<is_polyhedron<EncompassingType>::value, bool>
polytopeIsValid(const EncompassingType* a_polytope) {
  return true;
}

template <class EncompassingType>
enable_if_t<is_polygon<EncompassingType>::value, bool> polytopeIsValid(
    const EncompassingType* a_polytope) {
  return a_polytope->calculateAbsoluteVolume() <
             global_constants::MINIMUM_SURFACE_AREA_TO_TRACK ||
         a_polytope->getPlaneOfExistence().normal().calculateMagnitude() >
             0.999;
}

}  // namespace generic_cutting_details

//******************************************************************* //
//     Getting a polyhedrons volume or normalized moments
//******************************************************************* //
template <class ReturnType, class EncompassingType>
inline ReturnType getVolumeMoments(
    const EncompassingType& a_encompassing_polyhedron) {
  assert(generic_cutting_details::polytopeIsValid(a_encompassing_polyhedron));
  return ReturnType::calculateMoments(&a_encompassing_polyhedron);
}

template <class ReturnType, class CuttingMethod, class EncompassingType,
          class ReconstructionType>
__attribute__((pure)) __attribute__((hot)) inline ReturnType getVolumeMoments(
    const EncompassingType& a_encompassing_polyhedron,
    const ReconstructionType& a_reconstruction) {
  assert(generic_cutting_details::polytopeIsValid(a_encompassing_polyhedron));
  return generic_cutting_details::getVolumeMoments<
      ReturnType, CuttingMethod, EncompassingType, ReconstructionType>::
      getVolumeMomentsImplementation(a_encompassing_polyhedron,
                                     a_reconstruction);
}

template <class ReturnType, class CuttingMethod>
__attribute__((pure)) __attribute__((
    hot)) inline enable_if_t<IsOnlyVolume<ReturnType>::value, ReturnType>
getVolumeMoments(const Tet& a_encompassing_polyhedron,
                 const PlanarSeparator& a_reconstruction) {
  assert(generic_cutting_details::polytopeIsValid(a_encompassing_polyhedron));
  if (a_reconstruction.getNumberOfPlanes() == 1) {
    return getAnalyticVolume(a_encompassing_polyhedron, a_reconstruction[0]);
  } else {
    return generic_cutting_details::getVolumeMoments<ReturnType, CuttingMethod,
                                                     Tet, PlanarSeparator>::
        getVolumeMomentsImplementation(a_encompassing_polyhedron,
                                       a_reconstruction);
  }
}

template <class ReturnType, class CuttingMethod>
__attribute__((pure)) __attribute__((
    hot)) inline enable_if_t<IsOnlyVolume<ReturnType>::value, ReturnType>
getVolumeMoments(const RectangularCuboid& a_encompassing_polyhedron,
                 const PlanarSeparator& a_reconstruction) {
  assert(generic_cutting_details::polytopeIsValid(a_encompassing_polyhedron));
  if (a_reconstruction.getNumberOfPlanes() == 1) {
    return getAnalyticVolume(a_encompassing_polyhedron, a_reconstruction[0]);
  } else {
    return generic_cutting_details::getVolumeMoments<
        ReturnType, CuttingMethod, RectangularCuboid, PlanarSeparator>::
        getVolumeMomentsImplementation(a_encompassing_polyhedron,
                                       a_reconstruction);
  }
}

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ReconstructionType>
__attribute__((hot)) inline ReturnType getVolumeMoments(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction) {
  assert(generic_cutting_details::polytopeIsValid(a_polytope));
  return generic_cutting_details::getVolumeMomentsProvidedStorage<
      ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
      ReconstructionType>::getVolumeMomentsImplementation(a_polytope,
                                                          a_complete_polytope,
                                                          a_reconstruction);
}

//******************************************************************* //
//     Getting the normalized moments (or just volume) for a
//     polyhedron that may be localized and/or separated.
//******************************************************************* //
template <class ReturnType, class CuttingMethod, class EncompassingType,
          class ReconstructionType>
inline ReturnType getNormalizedVolumeMoments(
    const EncompassingType& a_encompassing_polyhedron,
    const ReconstructionType& a_reconstruction) {
  ReturnType object_to_return =
      IRL::getVolumeMoments<ReturnType, CuttingMethod>(
          a_encompassing_polyhedron, a_reconstruction);
  object_to_return.normalizeByVolume();
  return object_to_return;
}

//******************************************************************* //
//     Fraction of polyhedron's volume under reconstruction planes
//******************************************************************* //
template <class CuttingMethod, class EncompassingType, class ReconstructionType>
inline double getVolumeFraction(
    const EncompassingType& a_encompassing_polyhedron,
    const ReconstructionType& a_reconstruction) {
  double volume_fraction =
      getVolumeMoments<Volume, CuttingMethod>(a_encompassing_polyhedron,
                                              a_reconstruction) /
      safelyTiny(getVolumeMoments<Volume>(a_encompassing_polyhedron));
  return a_reconstruction.isFlipped() ? 1.0 - volume_fraction : volume_fraction;
}

namespace generic_cutting_details {

template <class ReturnType, class CuttingMethod, class EncompassingType>
ReturnType getVolumeMoments<
    ReturnType, CuttingMethod, EncompassingType, NullReconstruction,
    enable_if_t<IsNullReconstruction<NullReconstruction>::value>>::
    getVolumeMomentsImplementation(const EncompassingType& a_polytope,
                                   const NullReconstruction& a_reconstruction) {
  return ReturnType::calculateMoments(&a_polytope);
}

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType>
inline ReturnType getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    NullReconstruction,
    enable_if_t<IsNullReconstruction<NullReconstruction>::value>>::
    getVolumeMomentsImplementation(SegmentedPolytopeType* a_polytope,
                                   HalfEdgePolytopeType* a_complete_polytope,
                                   const NullReconstruction& a_reconstruction) {
  return ReturnType::calculateMoments(a_polytope);
}

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType>
inline ReturnType getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    Paraboloid, enable_if_t<IsParaboloidReconstruction<Paraboloid>::value>>::
    getVolumeMomentsImplementation(SegmentedPolytopeType* a_polytope,
                                   HalfEdgePolytopeType* a_complete_polytope,
                                   const Paraboloid& a_reconstruction) {
  return intersectPolyhedronWithParaboloid(a_polytope, a_complete_polytope,
                                           a_reconstruction);
}

template <class ReturnType, class CuttingMethod, class EncompassingType>
ReturnType getVolumeMoments<
    ReturnType, CuttingMethod, EncompassingType, PlanarSeparatorPathGroup,
    enable_if_t<IsPlanarSeparatorPathGroup<PlanarSeparatorPathGroup>::value>>::
    getVolumeMomentsImplementation(
        const EncompassingType& a_polytope,
        const PlanarSeparatorPathGroup& a_reconstruction) {
  return IRL::getVolumeMoments<ReturnType, CuttingMethod>(
      a_polytope, a_reconstruction.getFirstReconstruction());
}

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType>
inline ReturnType getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    PlanarSeparatorPathGroup,
    enable_if_t<IsPlanarSeparatorPathGroup<PlanarSeparatorPathGroup>::value>>::
    getVolumeMomentsImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const PlanarSeparatorPathGroup& a_reconstruction) {
  return IRL::getVolumeMoments<ReturnType, CuttingMethod>(
      a_polytope, a_complete_polytope,
      a_reconstruction.getFirstReconstruction());
}
//******************************************************************* //
//     Cut to obtain volume or VolumeMoments for volume under the planes.
//******************************************************************* //
template <class ReturnType, class CuttingMethod, class EncompassingType,
          class ReconstructionType>
ReturnType getVolumeMoments<
    ReturnType, CuttingMethod, EncompassingType, ReconstructionType,
    enable_if_t<isSimplexCutting<CuttingMethod>::value &&
                IsNotANullReconstruction<ReconstructionType>::value &&
                IsNotAPlanarSeparatorPathGroup<ReconstructionType>::value &&
                !(IsPlanarSeparator<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value)>>::
    getVolumeMomentsImplementation(const EncompassingType& a_polytope,
                                   const ReconstructionType& a_reconstruction) {
  return cutThroughSimplex<ReturnType>(a_polytope, a_reconstruction);
}

template <class ReturnType, class CuttingMethod, class EncompassingType,
          class ReconstructionType>
ReturnType getVolumeMoments<
    ReturnType, CuttingMethod, EncompassingType, ReconstructionType,
    enable_if_t<isHalfEdgeCutting<CuttingMethod>::value &&
                IsNotANullReconstruction<ReconstructionType>::value &&
                IsNotAPlanarSeparatorPathGroup<ReconstructionType>::value &&
                !(IsPlanarSeparator<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value)>>::
    getVolumeMomentsImplementation(const EncompassingType& a_polytope,
                                   const ReconstructionType& a_reconstruction) {
  return cutThroughHalfEdgeStructures<ReturnType>(a_polytope, a_reconstruction);
}

template <class ReturnType, class CuttingMethod, class EncompassingType,
          class ReconstructionType>
ReturnType getVolumeMoments<
    ReturnType, CuttingMethod, EncompassingType, ReconstructionType,
    enable_if_t<isRecursiveSimplexCutting<CuttingMethod>::value &&
                IsNotANullReconstruction<ReconstructionType>::value &&
                IsNotAPlanarSeparatorPathGroup<ReconstructionType>::value &&
                !(IsPlanarSeparator<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value)>>::
    getVolumeMomentsImplementation(const EncompassingType& a_polytope,
                                   const ReconstructionType& a_reconstruction) {
  return cutThroughRecursiveSimplex<ReturnType>(a_polytope, a_reconstruction);
}

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ReconstructionType>
inline ReturnType getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    ReconstructionType,
    enable_if_t<isHalfEdgeCutting<CuttingMethod>::value &&
                IsNotANullReconstruction<ReconstructionType>::value &&
                IsNotAPlanarSeparatorPathGroup<ReconstructionType>::value &&
                !(IsPlanarSeparator<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value)>>::
    getVolumeMomentsImplementation(SegmentedPolytopeType* a_polytope,
                                   HalfEdgePolytopeType* a_complete_polytope,
                                   const ReconstructionType& a_reconstruction) {
  ReturnType volume_moments;
  assert(a_polytope->checkValidHalfEdgeStructure());
  getVolumeMomentsForPolytope(a_polytope, a_complete_polytope, a_reconstruction,
                              &volume_moments);
  return volume_moments;
}

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ReconstructionType>
inline ReturnType getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    ReconstructionType,
    enable_if_t<isSimplexCutting<CuttingMethod>::value &&
                IsNotANullReconstruction<ReconstructionType>::value &&
                IsNotAPlanarSeparatorPathGroup<ReconstructionType>::value &&
                !(IsPlanarSeparator<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value)>>::
    getVolumeMomentsImplementation(SegmentedPolytopeType* a_polytope,
                                   HalfEdgePolytopeType* a_complete_polytope,
                                   const ReconstructionType& a_reconstruction) {
  ReturnType volume_moments;
  getVolumeMomentsForDecomposedPolytope(a_polytope, a_complete_polytope,
                                        a_reconstruction, &volume_moments);
  return volume_moments;
}

//******************************************************************* //
//     Getting the un-normalized SeparatedVolumeMomentsfor a
//     polyhedron that is separated.
//******************************************************************* //
template <class ReturnType, class CuttingMethod, class EncompassingType>
inline ReturnType
getVolumeMoments<ReturnType, CuttingMethod, EncompassingType, PlanarSeparator,
                 enable_if_t<is_separated_moments<ReturnType>::value>>::
    getVolumeMomentsImplementation(
        const EncompassingType& a_encompassing_polyhedron,
        const PlanarSeparator& a_separating_reconstruction) {
  typename ReturnType::moments_type volume_moments =
      IRL::getVolumeMoments<typename ReturnType::moments_type, CuttingMethod>(
          a_encompassing_polyhedron, a_separating_reconstruction);
  auto separated_volume_moments = ReturnType::fillWithComplementMoments(
      volume_moments, a_encompassing_polyhedron,
      a_separating_reconstruction.isFlipped());
  return separated_volume_moments;
}

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType>
inline ReturnType getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    PlanarSeparator,
    enable_if_t<isHalfEdgeCutting<CuttingMethod>::value &&
                is_separated_moments<ReturnType>::value>>::
    getVolumeMomentsImplementation(SegmentedPolytopeType* a_polytope,
                                   HalfEdgePolytopeType* a_complete_polytope,
                                   const PlanarSeparator& a_reconstruction) {
  typename ReturnType::moments_type encompassing_moments =
      ReturnType::moments_type::calculateMoments(a_polytope);
  auto volume_moments =
      IRL::getVolumeMoments<typename ReturnType::moments_type, HalfEdgeCutting>(
          a_polytope, a_complete_polytope, a_reconstruction);
  auto separated_volume_moments = ReturnType::fillWithComplementMoments(
      volume_moments, encompassing_moments, a_reconstruction.isFlipped());
  return separated_volume_moments;
}

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType>
inline ReturnType getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    PlanarSeparator,
    enable_if_t<isSimplexCutting<CuttingMethod>::value &&
                is_separated_moments<ReturnType>::value>>::
    getVolumeMomentsImplementation(SegmentedPolytopeType* a_polytope,
                                   HalfEdgePolytopeType* a_complete_polytope,
                                   const PlanarSeparator& a_reconstruction) {
  typename ReturnType::moments_type encompassing_moments =
      ReturnType::moments_type::calculateMoments(a_polytope);

  auto volume_moments =
      IRL::getVolumeMoments<typename ReturnType::moments_type, SimplexCutting>(
          a_polytope, a_complete_polytope, a_reconstruction);
  auto separated_volume_moments = ReturnType::fillWithComplementMoments(
      volume_moments, encompassing_moments, a_reconstruction.isFlipped());
  return separated_volume_moments;
}

}  // namespace generic_cutting_details

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_GENERIC_CUTTING_TPP_
