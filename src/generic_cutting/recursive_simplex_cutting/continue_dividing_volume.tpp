// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CONTINUE_DIVIDING_VOLUME_TPP_
#define SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CONTINUE_DIVIDING_VOLUME_TPP_

namespace IRL {

//******************************************************************* //
//******************************************************************* //

template <class SimplexType>
inline std::array<double, SimplexWrapper<SimplexType>::simplex_nvert>
getDistanceToVerticesFromPlane(const SimplexType& a_simplex,
                               const Plane& a_cutting_plane);

template <class ReconstructionType>
const Plane& getNextPlaneToCutBy(const ReconstructionType& a_reconstruction,
                                 const UnsignedIndex_t a_cutting_plane_index);

template <class VertexDistanceListType>
inline void flipDistanceToVertices(
    VertexDistanceListType* a_distance_to_vertices);

//******************************************************************* //
//     Function template definitions placed below this.               //
//******************************************************************* //
template <class SimplexType, class ReconstructionType>
std::array<double, SimplexWrapper<SimplexType>::simplex_nvert>
CutSimplexByNextPlaneAccessedByIndex::cutSimplexByNextPlane(
    const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
    const UnsignedIndex_t a_cutting_plane_index,
    LookupIndex_t* a_cutting_case_for_simplex_and_current_plane) {
  const auto& cutting_reconstruction =
      a_reconstruction.getCurrentReconstruction();

  const auto& cutting_plane =
      getNextPlaneToCutBy(cutting_reconstruction, a_cutting_plane_index);

  auto distances_to_simplex_vertices_from_cutting_plane =
      getDistanceToVerticesFromPlane(a_simplex, cutting_plane);
  if (cutting_reconstruction.isFlipped()) {
    flipDistanceToVertices(&distances_to_simplex_vertices_from_cutting_plane);
  }

  *a_cutting_case_for_simplex_and_current_plane =
      getGeometricCaseId(distances_to_simplex_vertices_from_cutting_plane);

  return distances_to_simplex_vertices_from_cutting_plane;
}

template <class SimplexType, class ReconstructionType>
inline UnsignedIndex_t findNextNonNullPlaneIndexToCutBy(
    const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
    const UnsignedIndex_t starting_plane_index) {
  static constexpr double small_number = DBL_EPSILON;
  const auto& cutting_reconstruction =
      a_reconstruction.getCurrentReconstruction();

  UnsignedIndex_t cutting_index_to_use = starting_plane_index;
  UnsignedIndex_t last_valid_index = cutting_reconstruction.getNumberOfPlanes();
  while (cutting_index_to_use != last_valid_index) {
    const auto& cutting_plane = cutting_reconstruction[cutting_index_to_use];
    for (const auto& pt : a_simplex) {
      if (cutting_reconstruction.flip() *
              cutting_plane.signedDistanceToPoint(pt) >
          small_number) {
        return cutting_index_to_use;
      }
    }
    ++cutting_index_to_use;
  }
  return cutting_index_to_use;  // Which will be last_valid_index if made it
                                // this far
}

template <class SimplexType, class ReconstructionType, class ReturnType>
bool earlyBranchBothFullyBelowAndFullyAboveSimplicesNoLink::
    earlyBranchingIfSimplexIsNotIntersectedByCurrentPlane(
        const SimplexType& a_simplex,
        const ReconstructionType& a_reconstruction,
        const UnsignedIndex_t a_cutting_plane_index,
        const LookupIndex_t a_cutting_case, ReturnType* a_moments_to_return) {
  if (SimplexWrapper<SimplexType>::isSimplexFullyBelowPlane(a_cutting_case)) {
    const UnsignedIndex_t next_plane_index_to_cut_by =
        findNextNonNullPlaneIndexToCutBy(a_simplex, a_reconstruction,
                                         a_cutting_plane_index + 1);
    getVolumeMomentsForSimplex(a_simplex, a_reconstruction,
                               next_plane_index_to_cut_by, a_moments_to_return);
    return true;
  }
  if (SimplexWrapper<SimplexType>::isSimplexFullyAbovePlane(a_cutting_case)) {
    return true;
  }
  return false;
}

template <class SimplexType, class ReconstructionType, class ReturnType>
bool earlyBranchBothFullyBelowAndFullyAboveSimplicesLink::
    earlyBranchingIfSimplexIsNotIntersectedByCurrentPlane(
        const SimplexType& a_simplex,
        const ReconstructionType& a_reconstruction,
        const UnsignedIndex_t a_cutting_plane_index,
        const LookupIndex_t a_cutting_case, ReturnType* a_moments_to_return) {
  if (SimplexWrapper<SimplexType>::isSimplexFullyBelowPlane(a_cutting_case)) {
    const UnsignedIndex_t next_plane_index_to_cut_by =
        findNextNonNullPlaneIndexToCutBy(a_simplex, a_reconstruction,
                                         a_cutting_plane_index + 1);
    getVolumeMomentsForSimplex(a_simplex, a_reconstruction,
                               next_plane_index_to_cut_by, a_moments_to_return);
    return true;
  }
  if (SimplexWrapper<SimplexType>::isSimplexFullyAbovePlane(a_cutting_case)) {
    if (a_reconstruction.hasNeighbor(a_cutting_plane_index)) {
      const auto& neighboring_reconstruction =
          a_reconstruction.getNeighbor(a_cutting_plane_index);
      getVolumeMomentsForSimplex(a_simplex, neighboring_reconstruction, 0,
                                 a_moments_to_return);
    }
    return true;
  }
  return false;
}

template <class SimplexType, class ReconstructionType, class VerticesList,
          class ReturnType>
void ContinueReducingVolumeToBeInternalToReconstruction::
    generateAndHandleSimplicesFromVolumeUnderPlane(
        const ReconstructionType& a_reconstruction,
        const LookupIndex_t a_cutting_case_for_simplex_and_current_plane,
        const VerticesList& a_simplex_vertices_and_intersection_points,
        const UnsignedIndex_t a_cutting_plane_index,
        ReturnType* a_moments_to_return) {
  for (LookupIndex_t t = 0;
       t < SimplexWrapper<SimplexType>::
               numberOfSimplicesInVolumeBelowPlaneAfterCutting(
                   a_cutting_case_for_simplex_and_current_plane);
       ++t) {
    const auto simplex_from_volume_below_plane =
        SimplexWrapper<SimplexType>::simplexFromCutSimplexVertices(
            a_simplex_vertices_and_intersection_points,
            a_cutting_case_for_simplex_and_current_plane, t);
    if (simplex_from_volume_below_plane.calculateAbsoluteVolume() >
        SimplexWrapper<SimplexType>::minimumAmountToTrack()) {
      getVolumeMomentsForSimplex(simplex_from_volume_below_plane,
                                 a_reconstruction, a_cutting_plane_index + 1,
                                 a_moments_to_return);
    }
  }
}

template <class SimplexType, class ReconstructionType, class VerticesList,
          class ReturnType>
void ShareVolumeAbovePlaneWithLinkedNeighbor::
    generateAndHandleSimplicesFromVolumeAbovePlane(
        const ReconstructionType& a_reconstruction,
        const LookupIndex_t a_cutting_case_for_simplex_and_current_plane,
        const VerticesList& a_simplex_vertices_and_intersection_points,
        const UnsignedIndex_t a_index_for_plane_that_just_cut,
        ReturnType* a_moments_to_return) {
  if (a_reconstruction.hasNeighbor(a_index_for_plane_that_just_cut)) {
    const auto& neighboring_reconstruction =
        a_reconstruction.getNeighbor(a_index_for_plane_that_just_cut);
    for (LookupIndex_t t = SimplexWrapper<SimplexType>::
             numberOfSimplicesInVolumeBelowPlaneAfterCutting(
                 a_cutting_case_for_simplex_and_current_plane);
         t < SimplexWrapper<SimplexType>::numberOfSimplicesInVolumeAfterCutting(
                 a_cutting_case_for_simplex_and_current_plane);
         ++t) {
      const auto simplex_from_volume_above_plane =
          SimplexWrapper<SimplexType>::simplexFromCutSimplexVertices(
              a_simplex_vertices_and_intersection_points,
              a_cutting_case_for_simplex_and_current_plane, t);
      if (simplex_from_volume_above_plane.calculateAbsoluteVolume() >
          SimplexWrapper<SimplexType>::minimumAmountToTrack()) {
        getVolumeMomentsForSimplex(simplex_from_volume_above_plane,
                                   neighboring_reconstruction, 0,
                                   a_moments_to_return);
      }
    }
  }
}

template <class ReconstructionType>
const Plane& getNextPlaneToCutBy(const ReconstructionType& a_reconstruction,
                                 const UnsignedIndex_t a_cutting_plane_index) {
  return a_reconstruction[a_cutting_plane_index];
}

template <class SimplexType>
std::array<double, SimplexWrapper<SimplexType>::simplex_nvert>
getDistanceToVerticesFromPlane(const SimplexType& a_simplex,
                               const Plane& a_cutting_plane) {
  std::array<double, SimplexWrapper<SimplexType>::simplex_nvert>
      distance_to_vertices;
  signedDistanceToVertices(a_simplex, a_cutting_plane, &distance_to_vertices);
  return distance_to_vertices;
}

template <class VertexDistanceListType>
void flipDistanceToVertices(VertexDistanceListType* a_distance_to_vertices) {
  for (auto& distance : (*a_distance_to_vertices)) {
    distance = -distance;
  }
}

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CONTINUE_DIVIDING_VOLUME_TPP_
