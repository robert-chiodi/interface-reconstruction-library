// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_GENERIC_CUTTING_H_
#define IRL_GENERIC_CUTTING_GENERIC_CUTTING_H_

#include <algorithm>
#include <type_traits>
#include <vector>

#include "irl/generic_cutting/general/class_classifications.h"
#include "irl/generic_cutting/generic_cutting_definitions.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_drivers.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_initializer.h"
#include "irl/generic_cutting/paraboloid_intersection/paraboloid_intersection.h"
#include "irl/generic_cutting/recursive_simplex_cutting/recursive_simplex_cutting_initializer.h"
#include "irl/generic_cutting/simplex_cutting/simplex_cutting_initializer.h"
#include "irl/helpers/SFINAE_boiler_plate.h"
#include "irl/moments/separated_volume_moments.h"
#include "irl/planar_reconstruction/null_reconstruction.h"
#include "irl/planar_reconstruction/planar_separator.h"
#include "irl/planar_reconstruction/planar_separator_path_group.h"

namespace IRL {

namespace generic_cutting_details {

template <class ReturnType, class CuttingMethod, class EncompassingType,
          class ReconstructionType, class Enable = void>
struct getVolumeMoments {};

template <class ReturnType, class CuttingMethod, class EncompassingType>
struct getVolumeMoments<
    ReturnType, CuttingMethod, EncompassingType, NullReconstruction,
    enable_if_t<IsNullReconstruction<NullReconstruction>::value>> {
  inline static ReturnType getVolumeMomentsImplementation(
      const EncompassingType& a_polytope,
      const NullReconstruction& a_reconstruction);
};

template <class ReturnType, class CuttingMethod, class EncompassingType>
struct getVolumeMoments<
    ReturnType, CuttingMethod, EncompassingType, PlanarSeparatorPathGroup,
    enable_if_t<IsPlanarSeparatorPathGroup<PlanarSeparatorPathGroup>::value>> {
  inline static ReturnType getVolumeMomentsImplementation(
      const EncompassingType& a_polytope,
      const PlanarSeparatorPathGroup& a_reconstruction);
};

template <class ReturnType, class CuttingMethod, class EncompassingType,
          class ReconstructionType>
struct getVolumeMoments<
    ReturnType, CuttingMethod, EncompassingType, ReconstructionType,
    enable_if_t<isSimplexCutting<CuttingMethod>::value &&
                IsNotANullReconstruction<ReconstructionType>::value &&
                IsNotAPlanarSeparatorPathGroup<ReconstructionType>::value &&
                !(IsPlanarSeparator<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value) &&
                !(IsParaboloidReconstruction<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value)>> {
  __attribute__((pure)) __attribute__((hot)) inline static ReturnType
  getVolumeMomentsImplementation(
      const EncompassingType& a_encompassing_polyhedron,
      const ReconstructionType& a_separating_reconstruction);
};

template <class ReturnType, class CuttingMethod, class EncompassingType,
          class ReconstructionType>
struct getVolumeMoments<
    ReturnType, CuttingMethod, EncompassingType, ReconstructionType,
    enable_if_t<isHalfEdgeCutting<CuttingMethod>::value &&
                IsNotANullReconstruction<ReconstructionType>::value &&
                IsNotAPlanarSeparatorPathGroup<ReconstructionType>::value &&
                !(IsPlanarSeparator<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value) &&
                !(IsParaboloidReconstruction<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value)>> {
  __attribute__((pure)) __attribute__((hot)) inline static ReturnType
  getVolumeMomentsImplementation(
      const EncompassingType& a_encompassing_polyhedron,
      const ReconstructionType& a_separating_reconstruction);
};

template <class ReturnType, class CuttingMethod, class EncompassingType,
          class ReconstructionType>
struct getVolumeMoments<
    ReturnType, CuttingMethod, EncompassingType, ReconstructionType,
    enable_if_t<isRecursiveSimplexCutting<CuttingMethod>::value &&
                IsNotANullReconstruction<ReconstructionType>::value &&
                IsNotAPlanarSeparatorPathGroup<ReconstructionType>::value &&
                !(IsPlanarSeparator<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value) &&
                !(IsParaboloidReconstruction<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value)>> {
  __attribute__((pure)) __attribute__((hot)) inline static ReturnType
  getVolumeMomentsImplementation(
      const EncompassingType& a_encompassing_polyhedron,
      const ReconstructionType& a_separating_reconstruction);
};

// Cut polyhedron for SeparatedMoments<VolumeMoments>
template <class ReturnType, class CuttingMethod, class EncompassingType>
struct getVolumeMoments<ReturnType, CuttingMethod, EncompassingType,
                        PlanarSeparator,
                        enable_if_t<is_separated_moments<ReturnType>::value>> {
  __attribute__((pure)) __attribute__((hot)) inline static ReturnType
  getVolumeMomentsImplementation(
      const EncompassingType& a_encompassing_polyhedron,
      const PlanarSeparator& a_separating_reconstruction);
};

template <class ReturnType, class CuttingMethod, class EncompassingType,
          class ScalarType>
struct getVolumeMoments<ReturnType, CuttingMethod, EncompassingType,
                        ParaboloidBase<ScalarType>,
                        enable_if_t<is_separated_moments<ReturnType>::value>> {
  __attribute__((pure)) __attribute__((hot)) inline static ReturnType
  getVolumeMomentsImplementation(
      const EncompassingType& a_encompassing_polyhedron,
      const ParaboloidBase<ScalarType>& a_separating_reconstruction);
};

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ReconstructionType,
          class Enable = void>
struct getVolumeMomentsProvidedStorage {};

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType>
struct getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    NullReconstruction,
    enable_if_t<IsNullReconstruction<NullReconstruction>::value>> {
  inline static ReturnType getVolumeMomentsImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const NullReconstruction& a_reconstruction);
};

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ScalarType>
struct getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    ParaboloidBase<ScalarType>,
    enable_if_t<IsParaboloidReconstruction<ParaboloidBase<ScalarType>>::value &&
                !is_separated_moments<ReturnType>::value>> {
  inline static ReturnType getVolumeMomentsImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ParaboloidBase<ScalarType>& a_reconstruction);
};

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType>
struct getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    PlanarSeparatorPathGroup,
    enable_if_t<IsPlanarSeparatorPathGroup<PlanarSeparatorPathGroup>::value>> {
  inline static ReturnType getVolumeMomentsImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const PlanarSeparatorPathGroup& a_reconstruction);
};

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ReconstructionType>
struct getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    ReconstructionType,
    enable_if_t<isHalfEdgeCutting<CuttingMethod>::value &&
                IsNotANullReconstruction<ReconstructionType>::value &&
                IsNotAPlanarSeparatorPathGroup<ReconstructionType>::value &&
                IsNotAParaboloidReconstruction<ReconstructionType>::value &&
                !(IsPlanarSeparator<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value) &&
                !(IsParaboloidReconstruction<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value)>> {
  __attribute__((hot)) inline static ReturnType getVolumeMomentsImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_reconstruction);
};

// Cut polyhedron for SeparatedMoments<VolumeMoments>
template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType>
struct getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    PlanarSeparator,
    enable_if_t<isHalfEdgeCutting<CuttingMethod>::value &&
                is_separated_moments<ReturnType>::value>> {
  __attribute__((hot)) inline static ReturnType getVolumeMomentsImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const PlanarSeparator& a_reconstruction);
};

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ScalarType>
struct getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    ParaboloidBase<ScalarType>,
    enable_if_t<isHalfEdgeCutting<CuttingMethod>::value &&
                is_separated_moments<ReturnType>::value>> {
  __attribute__((hot)) inline static ReturnType getVolumeMomentsImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ParaboloidBase<ScalarType>& a_reconstruction);
};

template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ReconstructionType>
struct getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    ReconstructionType,
    enable_if_t<isSimplexCutting<CuttingMethod>::value &&
                IsNotANullReconstruction<ReconstructionType>::value &&
                IsNotAPlanarSeparatorPathGroup<ReconstructionType>::value &&
                !(IsPlanarSeparator<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value) &&
                !(IsParaboloidReconstruction<ReconstructionType>::value &&
                  is_separated_moments<ReturnType>::value)>> {
  __attribute__((hot)) inline static ReturnType getVolumeMomentsImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_reconstruction);
};

// Cut polyhedron for SeparatedMoments<VolumeMoments>
template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType>
struct getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    PlanarSeparator,
    enable_if_t<isSimplexCutting<CuttingMethod>::value &&
                is_separated_moments<ReturnType>::value>> {
  __attribute__((hot)) inline static ReturnType getVolumeMomentsImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const PlanarSeparator& a_reconstruction);
};

// Cut polyhedron for SeparatedMoments<VolumeMoments>
template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ScalarType>
struct getVolumeMomentsProvidedStorage<
    ReturnType, CuttingMethod, SegmentedPolytopeType, HalfEdgePolytopeType,
    ParaboloidBase<ScalarType>,
    enable_if_t<isSimplexCutting<CuttingMethod>::value &&
                is_separated_moments<ReturnType>::value>> {
  __attribute__((hot)) inline static ReturnType getVolumeMomentsImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ParaboloidBase<ScalarType>& a_reconstruction);
};

}  // namespace generic_cutting_details

}  // namespace IRL

#include "irl/generic_cutting/generic_cutting.tpp"

#endif  // IRL_GENERIC_CUTTING_GENERIC_CUTTING_H_
