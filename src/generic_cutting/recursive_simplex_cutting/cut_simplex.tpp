// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CUT_SIMPLEX_TPP_
#define SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CUT_SIMPLEX_TPP_

namespace IRL {

namespace getVolumeMomentsForSimplexDetails {

template <class SimplexType, class ReconstructionType, class ReturnType,
          typename Enable = void>
struct getVolumeMomentsForSimplexStaticStructWrapper;

template <class SimplexType, class ReconstructionType, class ReturnType>
struct getVolumeMomentsForSimplexStaticStructWrapper<
    SimplexType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>> {
  static void getVolumeMomentsForSimplexImplementation(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      ReturnType* a_moments_to_return);
};

template <class SimplexType, class ReconstructionType, class ReturnType>
struct getVolumeMomentsForSimplexStaticStructWrapper<
    SimplexType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>> {
  static void getVolumeMomentsForSimplexImplementation(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      ReturnType* a_moments_to_return);
};

template <class SimplexType, class ReconstructionType, class ReturnType>
struct getVolumeMomentsForSimplexStaticStructWrapper<
    SimplexType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>> {
  static void getVolumeMomentsForSimplexImplementation(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      ReturnType* a_moments_to_return);
};

template <class SimplexType, class ReconstructionType, class ReturnType>
struct getVolumeMomentsForSimplexStaticStructWrapper<
    SimplexType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>> {
  static void getVolumeMomentsForSimplexImplementation(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      ReturnType* a_moments_to_return);
};

template <class SimplexType, class ReconstructionType, class ReturnType>
struct getVolumeMomentsForSimplexStaticStructWrapper<
    SimplexType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>> {
  static void getVolumeMomentsForSimplexImplementation(
      const SimplexType& a_simplex, const ReconstructionType& a_reconstruction,
      const UnsignedIndex_t a_cutting_plane_index,
      ReturnType* a_moments_to_return);
};

}  // namespace getVolumeMomentsForSimplexDetails

//******************************************************************* //
//     Function template definitions placed below this.
//******************************************************************* //
template <class SimplexType, class ReconstructionType, class ReturnType>
void getVolumeMomentsForSimplex(const SimplexType& a_simplex,
                                const ReconstructionType& a_reconstruction,
                                const UnsignedIndex_t a_cutting_plane_index,
                                ReturnType* a_moments_to_return) {
  getVolumeMomentsForSimplexDetails::
      getVolumeMomentsForSimplexStaticStructWrapper<
          SimplexType, ReconstructionType, ReturnType>::
          getVolumeMomentsForSimplexImplementation(a_simplex, a_reconstruction,
                                                   a_cutting_plane_index,
                                                   a_moments_to_return);
}

namespace getVolumeMomentsForSimplexDetails {

template <class SimplexType, class ReconstructionType, class ReturnType>
void getVolumeMomentsForSimplexStaticStructWrapper<
    SimplexType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>>::
    getVolumeMomentsForSimplexImplementation(
        const SimplexType& a_simplex,
        const ReconstructionType& a_reconstruction,
        const UnsignedIndex_t a_cutting_plane_index,
        ReturnType* a_moments_to_return) {
  DivideSimplexByPlanarReconstruction<
      AccumulateIntoScalar,
      KeepOnlyInternalReconstructionVolumeWithNoEarlyBranch>::
      execute(a_simplex, a_reconstruction, a_cutting_plane_index,
              a_moments_to_return);
}

template <class SimplexType, class ReconstructionType, class ReturnType>
void getVolumeMomentsForSimplexStaticStructWrapper<
    SimplexType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>>::
    getVolumeMomentsForSimplexImplementation(
        const SimplexType& a_simplex,
        const ReconstructionType& a_reconstruction,
        const UnsignedIndex_t a_cutting_plane_index,
        ReturnType* a_moments_to_return) {
  DivideSimplexByPlanarReconstruction<
      AccumulateIntoScalar,
      SpreadVolumeThroughLinkedLocalizerNetworkUsingEarlyBranch>::
      execute(a_simplex, a_reconstruction, a_cutting_plane_index,
              a_moments_to_return);
}

template <class SimplexType, class ReconstructionType, class ReturnType>
void getVolumeMomentsForSimplexStaticStructWrapper<
    SimplexType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>>::
    getVolumeMomentsForSimplexImplementation(
        const SimplexType& a_simplex,
        const ReconstructionType& a_reconstruction,
        const UnsignedIndex_t a_cutting_plane_index,
        ReturnType* a_moments_to_return) {
  DivideSimplexByPlanarReconstruction<
      PassToNestedType, KeepOnlyInternalReconstructionVolumeWithNoEarlyBranch>::
      execute(a_simplex, a_reconstruction, a_cutting_plane_index,
              a_moments_to_return);
}

template <class SimplexType, class ReconstructionType, class ReturnType>
void getVolumeMomentsForSimplexStaticStructWrapper<
    SimplexType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>>::
    getVolumeMomentsForSimplexImplementation(
        const SimplexType& a_simplex,
        const ReconstructionType& a_reconstruction,
        const UnsignedIndex_t a_cutting_plane_index,
        ReturnType* a_moments_to_return) {
  DivideSimplexByPlanarReconstruction<
      PassToNestedType,
      SpreadVolumeThroughLinkedLocalizerNetworkUsingEarlyBranch>::
      execute(a_simplex, a_reconstruction, a_cutting_plane_index,
              a_moments_to_return);
}

template <class SimplexType, class ReconstructionType, class ReturnType>
void getVolumeMomentsForSimplexStaticStructWrapper<
    SimplexType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>>::
    getVolumeMomentsForSimplexImplementation(
        const SimplexType& a_simplex,
        const ReconstructionType& a_reconstruction,
        const UnsignedIndex_t a_cutting_plane_index,
        ReturnType* a_moments_to_return) {
  DivideSimplexByPlanarReconstruction<
      AccumulateIntoCollection,
      SpreadVolumeThroughLinkedLocalizerNetworkUsingEarlyBranch>::
      execute(a_simplex, a_reconstruction, a_cutting_plane_index,
              a_moments_to_return);
}
}  // namespace getVolumeMomentsForSimplexDetails

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_RECURSIVE_SIMPLEX_CUTTING_CUT_SIMPLEX_TPP_
