// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_DRIVERS_TPP_
#define IRL_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_DRIVERS_TPP_

#include "irl/generic_cutting/general/encountered_pair_list.h"

namespace IRL {

// Foward declare getVolumeMoments to avoid circular dependency.
template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ReconstructionType>
__attribute__((hot)) inline ReturnType
getVolumeMoments(SegmentedPolytopeType *a_polytope,
                 HalfEdgePolytopeType *a_complete_polytope,
                 const ReconstructionType &a_reconstruction);

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
enable_if_t<HasAReconstructionLink<ReconstructionType>::value>
getVolumeMomentsForPolytope(SegmentedPolytopeType *a_polytope,
                            HalfEdgePolytopeType *a_complete_polytope,
                            const ReconstructionType &a_reconstruction,
                            EncounteredIdList *a_id_list,
                            ReturnType *a_moments_to_return);

namespace getVolumeMomentsForPolytopeDetails {
template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType, typename Enable = void>
struct getVolumeMomentsForPolytopeStaticStructWrapper;

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType *a_polytope,
      HalfEdgePolytopeType *a_complete_polytope,
      const ReconstructionType &a_reconstruction,
      ReturnType *a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType *a_polytope,
      HalfEdgePolytopeType *a_complete_polytope,
      const ReconstructionType &a_reconstruction,
      ReturnType *a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType *a_polytope,
      HalfEdgePolytopeType *a_complete_polytope,
      const ReconstructionType &a_reconstruction,
      ReturnType *a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType *a_polytope,
      HalfEdgePolytopeType *a_complete_polytope,
      const ReconstructionType &a_reconstruction,
      ReturnType *a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType *a_polytope,
      HalfEdgePolytopeType *a_complete_polytope,
      const ReconstructionType &a_reconstruction,
      ReturnType *a_moments_to_return);
};

// For ReconstructionLink cutting with informed origin.
template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType, typename Enable = void>
struct getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper;

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType *a_polytope,
      HalfEdgePolytopeType *a_complete_polytope,
      const ReconstructionType &a_reconstruction, EncounteredIdList *a_id_list,
      ReturnType *a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType *a_polytope,
      HalfEdgePolytopeType *a_complete_polytope,
      const ReconstructionType &a_reconstruction, EncounteredIdList *a_id_list,
      ReturnType *a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType *a_polytope,
      HalfEdgePolytopeType *a_complete_polytope,
      const ReconstructionType &a_reconstruction, EncounteredIdList *a_id_list,
      ReturnType *a_moments_to_return);
};

} // namespace getVolumeMomentsForPolytopeDetails

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType>
inline void
localizeInternalToReconstruction(SegmentedPolytopeType *a_polytope,
                                 HalfEdgePolytopeType *a_complete_polytope,
                                 const ReconstructionType &a_reconstruction);

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
inline void
splitAndShareThroughLinks(SegmentedPolytopeType *a_polytope,
                          HalfEdgePolytopeType *a_complete_polytope,
                          const ReconstructionType &a_reconstruction,
                          EncounteredIdList *a_id_list,
                          ReturnType *a_moments_to_return);

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void splitAndShareThroughLinks(SegmentedPolytopeType *a_polytope,
                               HalfEdgePolytopeType *a_complete_polytope,
                               const ReconstructionType &a_reconstruction,
                               EncounteredIdList *a_id_list,
                               ReturnType *a_moments_to_return);

//******************************************************************* //
//     Function template definitions placed below this
//******************************************************************* //

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytope(SegmentedPolytopeType *a_polytope,
                                 HalfEdgePolytopeType *a_complete_polytope,
                                 const ReconstructionType &a_reconstruction,
                                 ReturnType *a_moments_to_return) {
  getVolumeMomentsForPolytopeDetails::
      getVolumeMomentsForPolytopeStaticStructWrapper<
          SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType,
          ReturnType>::
          getVolumeMomentsForPolytopeImplementation(
              a_polytope, a_complete_polytope, a_reconstruction,
              a_moments_to_return);
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
enable_if_t<HasAReconstructionLink<ReconstructionType>::value>
getVolumeMomentsForPolytope(SegmentedPolytopeType *a_polytope,
                            HalfEdgePolytopeType *a_complete_polytope,
                            const ReconstructionType &a_reconstruction,
                            EncounteredIdList *a_id_list,
                            ReturnType *a_moments_to_return) {
  getVolumeMomentsForPolytopeDetails::
      getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
          SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType,
          ReturnType>::
          getVolumeMomentsForPolytopeImplementation(
              a_polytope, a_complete_polytope, a_reconstruction, a_id_list,
              a_moments_to_return);
}

namespace getVolumeMomentsForPolytopeDetails {

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType *a_polytope,
        HalfEdgePolytopeType *a_complete_polytope,
        const ReconstructionType &a_reconstruction,
        ReturnType *a_moments_to_return) {
  localizeInternalToReconstruction(a_polytope, a_complete_polytope,
                                   a_reconstruction);
  if (a_polytope->getNumberOfFaces() > 0) {
    *a_moments_to_return += IRL::getVolumeMoments<ReturnType, HalfEdgeCutting>(
        a_polytope, a_complete_polytope,
        a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType *a_polytope,
        HalfEdgePolytopeType *a_complete_polytope,
        const ReconstructionType &a_reconstruction,
        ReturnType *a_moments_to_return) {
  static EncounteredIdList id_list;
  splitAndShareThroughLinks(a_polytope, a_complete_polytope, a_reconstruction,
                            &id_list, a_moments_to_return);
  if (a_polytope->getNumberOfFaces() > 0) {
    *a_moments_to_return += IRL::getVolumeMoments<ReturnType, HalfEdgeCutting>(
        a_polytope, a_complete_polytope,
        a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType *a_polytope,
        HalfEdgePolytopeType *a_complete_polytope,
        const ReconstructionType &a_reconstruction,
        ReturnType *a_moments_to_return) {
  localizeInternalToReconstruction(a_polytope, a_complete_polytope,
                                   a_reconstruction);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  if (a_polytope->getNumberOfFaces() > 0) {
    *a_moments_to_return +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, HalfEdgeCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType *a_polytope,
        HalfEdgePolytopeType *a_complete_polytope,
        const ReconstructionType &a_reconstruction,
        ReturnType *a_moments_to_return) {
  static EncounteredIdList id_list;
  splitAndShareThroughLinks(a_polytope, a_complete_polytope, a_reconstruction,
                            &id_list, a_moments_to_return);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;

  if (a_polytope->getNumberOfFaces() > 0) {
    *a_moments_to_return +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, HalfEdgeCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType *a_polytope,
        HalfEdgePolytopeType *a_complete_polytope,
        const ReconstructionType &a_reconstruction,
        ReturnType *a_moments_to_return) {
  static EncounteredIdList id_list;
  splitAndShareThroughLinks(a_polytope, a_complete_polytope, a_reconstruction,
                            &id_list, a_moments_to_return);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;

  if (a_polytope->getNumberOfFaces() > 0) {
    assert(a_reconstruction.isIdSet());
    (*a_moments_to_return)[a_reconstruction.getId()] +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, HalfEdgeCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType *a_polytope,
        HalfEdgePolytopeType *a_complete_polytope,
        const ReconstructionType &a_reconstruction,
        EncounteredIdList *a_id_list, ReturnType *a_moments_to_return) {
  splitAndShareThroughLinks(a_polytope, a_complete_polytope, a_reconstruction,
                            a_id_list, a_moments_to_return);
  if (a_polytope->getNumberOfFaces() > 0) {
    *a_moments_to_return += IRL::getVolumeMoments<ReturnType, HalfEdgeCutting>(
        a_polytope, a_complete_polytope,
        a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType *a_polytope,
        HalfEdgePolytopeType *a_complete_polytope,
        const ReconstructionType &a_reconstruction,
        EncounteredIdList *a_id_list, ReturnType *a_moments_to_return) {
  splitAndShareThroughLinks(a_polytope, a_complete_polytope, a_reconstruction,
                            a_id_list, a_moments_to_return);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  if (a_polytope->getNumberOfFaces() > 0) {
    *a_moments_to_return +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, HalfEdgeCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType *a_polytope,
        HalfEdgePolytopeType *a_complete_polytope,
        const ReconstructionType &a_reconstruction,
        EncounteredIdList *a_id_list, ReturnType *a_moments_to_return) {
  splitAndShareThroughLinks(a_polytope, a_complete_polytope, a_reconstruction,
                            a_id_list, a_moments_to_return);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  if (a_polytope->getNumberOfFaces() > 0) {
    assert(a_reconstruction.isIdSet());
    (*a_moments_to_return)[a_reconstruction.getId()] +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, HalfEdgeCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

} // namespace getVolumeMomentsForPolytopeDetails

template <class SegmentedHalfEdgePolytopeType>
static constexpr enable_if_t<
    is_polyhedron<SegmentedHalfEdgePolytopeType>::value, double>
getGlobalVolumeTolerance(void) {
  return global_constants::MINIMUM_VOLUME_TO_TRACK;
}

template <class SegmentedHalfEdgePolytopeType>
static constexpr enable_if_t<is_polygon<SegmentedHalfEdgePolytopeType>::value,
                             double>
getGlobalVolumeTolerance(void) {
  return global_constants::MINIMUM_SURFACE_AREA_TO_TRACK;
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType>
void localizeInternalToReconstruction(
    SegmentedPolytopeType *a_polytope,
    HalfEdgePolytopeType *a_complete_polytope,
    const ReconstructionType &a_reconstruction) {
  const auto &cutting_reconstruction =
      a_reconstruction.getCurrentReconstruction();

  for (const auto &plane : cutting_reconstruction) {
    const auto cutting_plane = cutting_reconstruction.isFlipped()
                                   ? plane.generateFlippedPlane()
                                   : plane;
    if (a_polytope->getNumberOfFaces() == 0) {
      return;
    }
    truncateHalfEdgePolytope(a_polytope, a_complete_polytope, cutting_plane);
  }
}

namespace details {

// Check if ID is present but take advantage of the fact we know the last id
// is the cell's own id (and is therefore ignored)
inline bool isIdPresent(const EncounteredIdList &a_list,
                        const UnsignedIndex_t a_id) {
  const auto end = a_list.end() - 1;
  return std::find(a_list.begin(), end, a_id) != end;
}

} // namespace details

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void splitAndShareThroughLinks(SegmentedPolytopeType *a_polytope,
                               HalfEdgePolytopeType *a_complete_polytope,
                               const ReconstructionType &a_reconstruction,
                               EncounteredIdList *a_id_list,
                               ReturnType *a_moments_to_return) {
  const auto &cutting_reconstruction =
      a_reconstruction.getCurrentReconstruction();
  a_id_list->push_back(a_reconstruction.getId());

  SegmentedPolytopeType clipped_polytope;
  for (UnsignedIndex_t plane_index = 0;
       plane_index < cutting_reconstruction.getNumberOfPlanes();
       ++plane_index) {
    if (a_polytope->getNumberOfFaces() == 0) {
      break;
    }

    const auto cutting_plane =
        cutting_reconstruction.isFlipped()
            ? cutting_reconstruction[plane_index].generateFlippedPlane()
            : cutting_reconstruction[plane_index];

    if (a_reconstruction.hasNeighbor(plane_index)) {
      if (!details::isIdPresent(
              *a_id_list, a_reconstruction.getNeighbor(plane_index).getId())) {
        splitHalfEdgePolytope(a_polytope, &clipped_polytope,
                              a_complete_polytope, cutting_plane);
        if (clipped_polytope.getNumberOfFaces() > 0) {
          getVolumeMomentsForPolytope(&clipped_polytope, a_complete_polytope,
                                      a_reconstruction.getNeighbor(plane_index),
                                      a_id_list, a_moments_to_return);
        }
      }
    } else {
      truncateHalfEdgePolytope(a_polytope, a_complete_polytope, cutting_plane);
    }
  }
  a_id_list->pop_back();
}

} // namespace IRL

#endif // SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_DRIVERS_TPP_
