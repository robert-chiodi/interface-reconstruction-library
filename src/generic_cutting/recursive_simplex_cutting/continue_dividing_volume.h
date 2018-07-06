// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CONTINUE_DIVIDING_VOLUME_H_
#define SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CONTINUE_DIVIDING_VOLUME_H_

#include "src/generic_cutting/general/class_classifications.h"
#include "src/generic_cutting/recursive_simplex_cutting/simplex_wrapper.h"
#include "src/helpers/SFINAE_boiler_plate.h"
#include "src/parameters/defined_types.h"

namespace IRL {
// Forward declaration
template <class SimplexType, class ReconstructionType, class ReturnType>
void getVolumeMomentsForSimplex(const SimplexType& a_simplex,
                                const ReconstructionType& a_reconstruction,
                                const UnsignedIndex_t a_cutting_plane_index,
                                ReturnType* a_moments_to_return);

// This is the main structure used in cut_simplex_building_blocks.
// Can select different functions for HowToCutSimplexByNextPlane,
// HowToGenerateAndHandleSimplicesFromVolumeUnderPlane, and
// HowToGenerateAndHandleSimplicesFromVolumeAbovePlane, which are
// structs that wrap static functions .
template <class HowToCutSimplexByNextPlane,
          class HowToBranchEarlyIfNoIntersection,
          class HowToGenerateAndHandleSimplicesFromVolumeUnderPlane,
          class HowToGenerateAndHandleSimplicesFromVolumeAbovePlane>
class ContinueDividingVolumeByPlanarReconstruction {
 public:
  template <class SimplexType, class ReconstructionType, class ReturnType>
  static void continueDividingVolumeByPlanarReconstruction(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      ReturnType* a_moments_to_return) {
    assert(a_cutting_plane_index <
           a_reconstruction.getCurrentReconstruction().getNumberOfPlanes());
    LookupIndex_t cutting_case_for_simplex_and_current_plane;
    const auto distances_to_simplex_vertices_from_cutting_plane =
        cutSimplexByNextPlane(a_simplex, a_reconstruction,
                              a_cutting_plane_index,
                              &cutting_case_for_simplex_and_current_plane);

    const bool took_early_branch =
        earlyBranchingIfSimplexIsNotIntersectedByCurrentPlane(
            a_simplex, a_reconstruction, a_cutting_plane_index,
            cutting_case_for_simplex_and_current_plane, a_moments_to_return);
    if (!took_early_branch) {
      const auto simplex_vertices_and_intersection_points =
          SimplexWrapper<SimplexType>::getSimplexVerticesAndIntersectionPoints(
              a_simplex, distances_to_simplex_vertices_from_cutting_plane,
              cutting_case_for_simplex_and_current_plane);
      generateAndHandleSimplicesFromVolumeUnderPlane<SimplexType>(
          a_reconstruction, cutting_case_for_simplex_and_current_plane,
          simplex_vertices_and_intersection_points, a_cutting_plane_index,
          a_moments_to_return);
      generateAndHandleSimplicesFromVolumeAbovePlane<SimplexType>(
          a_reconstruction, cutting_case_for_simplex_and_current_plane,
          simplex_vertices_and_intersection_points, a_cutting_plane_index,
          a_moments_to_return);
    }
  }

 private:
  template <class SimplexType, class ReconstructionType>
  static std::array<double, SimplexWrapper<SimplexType>::simplex_nvert>
  cutSimplexByNextPlane(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      LookupIndex_t* a_cutting_case_for_simplex_and_current_plane) {
    return HowToCutSimplexByNextPlane::cutSimplexByNextPlane(
        a_simplex, a_reconstruction, a_cutting_plane_index,
        a_cutting_case_for_simplex_and_current_plane);
  }

  template <class SimplexType, class ReconstructionType, class ReturnType>
  static bool earlyBranchingIfSimplexIsNotIntersectedByCurrentPlane(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      const LookupIndex_t a_cutting_case, ReturnType* a_moments_to_return) {
    return HowToBranchEarlyIfNoIntersection::
        earlyBranchingIfSimplexIsNotIntersectedByCurrentPlane(
            a_simplex, a_reconstruction, a_cutting_plane_index, a_cutting_case,
            a_moments_to_return);
  }

  template <class SimplexType, class ReconstructionType, class VerticesList,
            class ReturnType>
  static void generateAndHandleSimplicesFromVolumeUnderPlane(
      const ReconstructionType& a_reconstruction,
      const LookupIndex_t a_cutting_case_for_simplex_and_current_plane,
      const VerticesList& a_simplex_vertices_and_intersection_points,
      const UnsignedIndex_t a_cutting_plane_index,
      ReturnType* a_moments_to_return) {
    HowToGenerateAndHandleSimplicesFromVolumeUnderPlane::
        template generateAndHandleSimplicesFromVolumeUnderPlane<SimplexType>(
            a_reconstruction, a_cutting_case_for_simplex_and_current_plane,
            a_simplex_vertices_and_intersection_points, a_cutting_plane_index,
            a_moments_to_return);
  }

  template <class SimplexType, class ReconstructionType, class VerticesList,
            class ReturnType>
  static void generateAndHandleSimplicesFromVolumeAbovePlane(
      const ReconstructionType& a_reconstruction,
      const LookupIndex_t a_cutting_case_for_simplex_and_current_plane,
      const VerticesList& a_simplex_vertices_and_intersection_points,
      const UnsignedIndex_t a_index_for_plane_that_just_cut,
      ReturnType* a_moments_to_return) {
    HowToGenerateAndHandleSimplicesFromVolumeAbovePlane::
        template generateAndHandleSimplicesFromVolumeAbovePlane<SimplexType>(
            a_reconstruction, a_cutting_case_for_simplex_and_current_plane,
            a_simplex_vertices_and_intersection_points,
            a_index_for_plane_that_just_cut, a_moments_to_return);
  }
};

struct CutSimplexByNextPlaneAccessedByIndex {
  template <class SimplexType, class ReconstructionType>
  inline static std::array<double, SimplexWrapper<SimplexType>::simplex_nvert>
  cutSimplexByNextPlane(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      LookupIndex_t* a_cutting_case_for_simplex_and_current_plane);
};

struct ContinueReducingVolumeToBeInternalToReconstruction {
  template <class SimplexType, class ReconstructionType, class VerticesList,
            class ReturnType>
  inline static void generateAndHandleSimplicesFromVolumeUnderPlane(
      const ReconstructionType& a_reconstruction,
      const LookupIndex_t a_cutting_case_for_simplex_and_current_plane,
      const VerticesList& a_simplex_vertices_and_intersection_points,
      const UnsignedIndex_t a_cutting_plane_index,
      ReturnType* a_moments_to_return);
};

struct doNoEarlyBranching {
  template <class SimplexType, class ReconstructionType, class ReturnType>
  inline static constexpr bool
  earlyBranchingIfSimplexIsNotIntersectedByCurrentPlane(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      const LookupIndex_t a_cutting_case, ReturnType* a_moments_to_return) {
    return false;
  }
};

// Handle cases with reconstruction links different from those without
struct earlyBranchBothFullyBelowAndFullyAboveSimplicesNoLink {
  template <class SimplexType, class ReconstructionType, class ReturnType>
  inline static bool earlyBranchingIfSimplexIsNotIntersectedByCurrentPlane(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      const LookupIndex_t a_cutting_case, ReturnType* a_moments_to_return);
};

struct earlyBranchBothFullyBelowAndFullyAboveSimplicesLink {
  template <class SimplexType, class ReconstructionType, class ReturnType>
  inline static bool earlyBranchingIfSimplexIsNotIntersectedByCurrentPlane(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      const LookupIndex_t a_cutting_case, ReturnType* a_moments_to_return);
};

struct IgnoreVolumeAbovePlane {
  template <class SimplexType, class ReconstructionType, class VerticesList,
            class ReturnType>
  inline static void generateAndHandleSimplicesFromVolumeAbovePlane(
      const ReconstructionType& a_reconstruction,
      const LookupIndex_t a_cutting_case_for_simplex_and_current_plane,
      const VerticesList& a_simplex_vertices_and_intersection_points,
      const UnsignedIndex_t a_index_for_plane_that_just_cut,
      ReturnType* a_moments_to_return) {}
};

struct ShareVolumeAbovePlaneWithLinkedNeighbor {
  template <class SimplexType, class ReconstructionType, class VerticesList,
            class ReturnType>
  inline static void generateAndHandleSimplicesFromVolumeAbovePlane(
      const ReconstructionType& a_reconstruction,
      const LookupIndex_t a_cutting_case_for_simplex_and_current_plane,
      const VerticesList& a_simplex_vertices_and_intersection_points,
      const UnsignedIndex_t a_index_for_plane_that_just_cut,
      ReturnType* a_moments_to_return);
};

//******************************************************************* //
//     Simplifying definitions of common treatments of cut volumes.   //
//******************************************************************* //
using KeepOnlyInternalReconstructionVolumeWithNoEarlyBranch =
    ContinueDividingVolumeByPlanarReconstruction<
        CutSimplexByNextPlaneAccessedByIndex, doNoEarlyBranching,
        ContinueReducingVolumeToBeInternalToReconstruction,
        IgnoreVolumeAbovePlane>;

using SpreadVolumeThroughLinkedLocalizerNetworkWithNoEarlyBranch =
    ContinueDividingVolumeByPlanarReconstruction<
        CutSimplexByNextPlaneAccessedByIndex, doNoEarlyBranching,
        ContinueReducingVolumeToBeInternalToReconstruction,
        ShareVolumeAbovePlaneWithLinkedNeighbor>;

using KeepOnlyInternalReconstructionVolumeUsingEarlyBranch =
    ContinueDividingVolumeByPlanarReconstruction<
        CutSimplexByNextPlaneAccessedByIndex,
        earlyBranchBothFullyBelowAndFullyAboveSimplicesNoLink,
        ContinueReducingVolumeToBeInternalToReconstruction,
        IgnoreVolumeAbovePlane>;

using SpreadVolumeThroughLinkedLocalizerNetworkUsingEarlyBranch =
    ContinueDividingVolumeByPlanarReconstruction<
        CutSimplexByNextPlaneAccessedByIndex,
        earlyBranchBothFullyBelowAndFullyAboveSimplicesLink,
        ContinueReducingVolumeToBeInternalToReconstruction,
        ShareVolumeAbovePlaneWithLinkedNeighbor>;

}  // namespace IRL

#include "src/generic_cutting/recursive_simplex_cutting/continue_dividing_volume.tpp"

#endif  // SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CONTINUE_DIVIDING_VOLUME_H_
