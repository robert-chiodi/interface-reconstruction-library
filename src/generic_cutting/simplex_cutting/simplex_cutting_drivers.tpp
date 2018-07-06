// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_SIMPLEX_CUTTING_SIMPLEX_CUTTING_DRIVERS_TPP_
#define SRC_GENERIC_CUTTING_SIMPLEX_CUTTING_SIMPLEX_CUTTING_DRIVERS_TPP_

#include "src/generic_cutting/simplex_cutting/simplex_cutting.h"

namespace IRL {

// Foward declare getVolumeMoments to avoid circular dependency.
template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ReconstructionType>
__attribute__((hot)) inline ReturnType getVolumeMoments(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction);

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
enable_if_t<HasAReconstructionLink<ReconstructionType>::value>
getVolumeMomentsForDecomposedPolytope(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_originating_reconstruction,
    const ReconstructionType& a_reconstruction,
    ReturnType* a_moments_to_return);

namespace getVolumeMomentsForDecomposedPolytopeDetails {
template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType, typename Enable = void>
struct getVolumeMomentsForDecomposedPolytopeStaticStructWrapper;

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForDecomposedPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>> {
  static void getVolumeMomentsForDecomposedPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForDecomposedPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>> {
  static void getVolumeMomentsForDecomposedPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForDecomposedPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>> {
  static void getVolumeMomentsForDecomposedPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForDecomposedPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>> {
  static void getVolumeMomentsForDecomposedPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForDecomposedPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>> {
  static void getVolumeMomentsForDecomposedPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

// For ReconstructionLink cutting with informed origin.
template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType, typename Enable = void>
struct getVolumeMomentsForDecomposedPolytopeKnownOriginStaticStructWrapper;

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForDecomposedPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>> {
  static void getVolumeMomentsForDecomposedPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_originating_reconstruction,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForDecomposedPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>> {
  static void getVolumeMomentsForDecomposedPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_originating_reconstruction,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForDecomposedPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>> {
  static void getVolumeMomentsForDecomposedPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_originating_reconstruction,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

}  // namespace getVolumeMomentsForDecomposedPolytopeDetails

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType>
inline void localizeSimplexInternalToReconstruction(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction);

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
inline void splitAndShareSimplexThroughLinks(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction,
    ReturnType* a_moments_to_return);

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void splitAndShareSimplexThroughLinks(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_originating_reconstruction,
    const ReconstructionType& a_reconstruction,
    ReturnType* a_moments_to_return);

//******************************************************************* //
//     Function template definitions placed below this
//******************************************************************* //

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForDecomposedPolytope(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction,
    ReturnType* a_moments_to_return) {
  getVolumeMomentsForDecomposedPolytopeDetails::
      getVolumeMomentsForDecomposedPolytopeStaticStructWrapper<
          SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType,
          ReturnType>::
          getVolumeMomentsForDecomposedPolytopeImplementation(
              a_polytope, a_complete_polytope, a_reconstruction,
              a_moments_to_return);
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
enable_if_t<HasAReconstructionLink<ReconstructionType>::value>
getVolumeMomentsForDecomposedPolytope(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_originating_reconstruction,
    const ReconstructionType& a_reconstruction,
    ReturnType* a_moments_to_return) {
  getVolumeMomentsForDecomposedPolytopeDetails::
      getVolumeMomentsForDecomposedPolytopeKnownOriginStaticStructWrapper<
          SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType,
          ReturnType>::
          getVolumeMomentsForDecomposedPolytopeImplementation(
              a_polytope, a_complete_polytope, a_originating_reconstruction,
              a_reconstruction, a_moments_to_return);
}

namespace getVolumeMomentsForDecomposedPolytopeDetails {

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForDecomposedPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>>::
    getVolumeMomentsForDecomposedPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  localizeSimplexInternalToReconstruction(a_polytope, a_complete_polytope,
                                          a_reconstruction);
  if (a_polytope->getNumberOfSimplicesInDecomposition() > 0) {
    *a_moments_to_return += IRL::getVolumeMoments<ReturnType, SimplexCutting>(
        a_polytope, a_complete_polytope,
        a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForDecomposedPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>>::
    getVolumeMomentsForDecomposedPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  splitAndShareSimplexThroughLinks(a_polytope, a_complete_polytope,
                                   a_reconstruction, a_moments_to_return);
  if (a_polytope->getNumberOfSimplicesInDecomposition() > 0) {
    *a_moments_to_return += IRL::getVolumeMoments<ReturnType, SimplexCutting>(
        a_polytope, a_complete_polytope,
        a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForDecomposedPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>>::
    getVolumeMomentsForDecomposedPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  localizeSimplexInternalToReconstruction(a_polytope, a_complete_polytope,
                                          a_reconstruction);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  if (a_polytope->getNumberOfSimplicesInDecomposition() > 0) {
    *a_moments_to_return +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, SimplexCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForDecomposedPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>>::
    getVolumeMomentsForDecomposedPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  splitAndShareSimplexThroughLinks(a_polytope, a_complete_polytope,
                                   a_reconstruction, a_moments_to_return);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  if (a_polytope->getNumberOfSimplicesInDecomposition() > 0) {
    *a_moments_to_return +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, SimplexCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForDecomposedPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>>::
    getVolumeMomentsForDecomposedPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  splitAndShareSimplexThroughLinks(a_polytope, a_complete_polytope,
                                   a_reconstruction, a_moments_to_return);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  if (a_polytope->getNumberOfSimplicesInDecomposition() > 0) {
    assert(a_reconstruction.isIdSet());
    (*a_moments_to_return)[a_reconstruction.getId()] +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, SimplexCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForDecomposedPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>>::
    getVolumeMomentsForDecomposedPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_originating_reconstruction,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  splitAndShareSimplexThroughLinks(a_polytope, a_complete_polytope,
                                   a_originating_reconstruction,
                                   a_reconstruction, a_moments_to_return);
  if (a_polytope->getNumberOfSimplicesInDecomposition() > 0) {
    *a_moments_to_return += IRL::getVolumeMoments<ReturnType, SimplexCutting>(
        a_polytope, a_complete_polytope,
        a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForDecomposedPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>>::
    getVolumeMomentsForDecomposedPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_originating_reconstruction,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  splitAndShareSimplexThroughLinks(a_polytope, a_complete_polytope,
                                   a_originating_reconstruction,
                                   a_reconstruction, a_moments_to_return);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  if (a_polytope->getNumberOfSimplicesInDecomposition() > 0) {
    *a_moments_to_return +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, SimplexCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForDecomposedPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>>::
    getVolumeMomentsForDecomposedPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_originating_reconstruction,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  splitAndShareSimplexThroughLinks(a_polytope, a_complete_polytope,
                                   a_originating_reconstruction,
                                   a_reconstruction, a_moments_to_return);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  if (a_polytope->getNumberOfSimplicesInDecomposition() > 0) {
    assert(a_reconstruction.isIdSet());
    (*a_moments_to_return)[a_reconstruction.getId()] +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, SimplexCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

}  // namespace getVolumeMomentsForDecomposedPolytopeDetails

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType>
void localizeSimplexInternalToReconstruction(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction) {
  const auto& cutting_reconstruction =
      a_reconstruction.getCurrentReconstruction();

  for (const auto& plane : cutting_reconstruction) {
    const auto cutting_plane = cutting_reconstruction.isFlipped()
                                   ? plane.generateFlippedPlane()
                                   : plane;

    if (a_polytope->getNumberOfSimplicesInDecomposition() == 0) {
      break;
    }
    truncateDecomposedPolytope(a_polytope, a_complete_polytope, cutting_plane);
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void splitAndShareSimplexThroughLinks(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction,
    ReturnType* a_moments_to_return) {
  const auto& cutting_reconstruction =
      a_reconstruction.getCurrentReconstruction();

  SegmentedPolytopeType clipped_polytope(*a_complete_polytope);
  for (UnsignedIndex_t plane_index = 0;
       plane_index < cutting_reconstruction.getNumberOfPlanes();
       ++plane_index) {
    const auto cutting_plane =
        cutting_reconstruction.isFlipped()
            ? cutting_reconstruction[plane_index].generateFlippedPlane()
            : cutting_reconstruction[plane_index];
    if (a_polytope->getNumberOfSimplicesInDecomposition() == 0) {
      break;
    }

    if (a_reconstruction.hasNeighbor(plane_index)) {
      clipped_polytope.clear();
      splitDecomposedPolytope(a_polytope, &clipped_polytope,
                              a_complete_polytope, cutting_plane);
      getVolumeMomentsForDecomposedPolytope(
          &clipped_polytope, a_complete_polytope, a_reconstruction,
          a_reconstruction.getNeighbor(plane_index), a_moments_to_return);
    } else {
      truncateDecomposedPolytope(a_polytope, a_complete_polytope,
                                 cutting_plane);
    }
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void splitAndShareSimplexThroughLinks(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_originating_reconstruction,
    const ReconstructionType& a_reconstruction,
    ReturnType* a_moments_to_return) {
  const auto& cutting_reconstruction =
      a_reconstruction.getCurrentReconstruction();

  SegmentedPolytopeType clipped_polytope(*a_complete_polytope);
  for (UnsignedIndex_t plane_index = 0;
       plane_index < cutting_reconstruction.getNumberOfPlanes();
       ++plane_index) {
    if (a_reconstruction.getNeighborAddress(plane_index) ==
        a_originating_reconstruction.getLinkingReconstructionAddress()) {
      continue;
    }
    const auto cutting_plane =
        cutting_reconstruction.isFlipped()
            ? cutting_reconstruction[plane_index].generateFlippedPlane()
            : cutting_reconstruction[plane_index];
    if (a_polytope->getNumberOfSimplicesInDecomposition() == 0) {
      break;
    }

    if (a_reconstruction.hasNeighbor(plane_index)) {
      clipped_polytope.clear();
      splitDecomposedPolytope(a_polytope, &clipped_polytope,
                              a_complete_polytope, cutting_plane);
      getVolumeMomentsForDecomposedPolytope(
          &clipped_polytope, a_complete_polytope, a_reconstruction,
          a_reconstruction.getNeighbor(plane_index), a_moments_to_return);
    } else {
      truncateDecomposedPolytope(a_polytope, a_complete_polytope,
                                 cutting_plane);
    }
  }
}

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_SIMPLEX_CUTTING_SIMPLEX_CUTTING_DRIVERS_TPP_
