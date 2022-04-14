// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_GENERAL_CLASS_CLASSIFICATIONS_H_
#define IRL_GENERIC_CUTTING_GENERAL_CLASS_CLASSIFICATIONS_H_

#include "irl/generic_cutting/default_cutting_method.h"
#include "irl/moments/moments_type_traits.h"
#include "irl/paraboloid_reconstruction/paraboloid_reconstruction_type_traits.h"
#include "irl/planar_reconstruction/planar_reconstruction_type_traits.h"

namespace IRL {

// Cutting methods
template <class C>
struct isHalfEdgeCutting : std::false_type {};

template <class C>
struct isHalfEdgeCutting<const C> : isHalfEdgeCutting<C> {};

template <>
struct isHalfEdgeCutting<HalfEdgeCutting> : std::true_type {};

template <class C>
struct isRecursiveSimplexCutting : std::false_type {};

template <class C>
struct isRecursiveSimplexCutting<const C> : isRecursiveSimplexCutting<C> {};

template <>
struct isRecursiveSimplexCutting<RecursiveSimplexCutting> : std::true_type {};

template <class C>
struct isSimplexCutting : std::false_type {};

template <class C>
struct isSimplexCutting<const C> : isSimplexCutting<C> {};

template <>
struct isSimplexCutting<SimplexCutting> : std::true_type {};

// More explanatory characterizations of different classes used for cutting.
template <class ReconstructionType>
struct IsNullReconstruction {
  static constexpr bool value =
      has_null_reconstruction<ReconstructionType>::value;
};

template <class ReconstructionType>
struct IsNotANullReconstruction {
  static constexpr bool value =
      !IsNullReconstruction<ReconstructionType>::value;
};

template <class ReconstructionType>
struct HasAParaboloidReconstruction {
  static constexpr bool value =
      has_paraboloid_reconstruction<ReconstructionType>::value;
};

template <class ReconstructionType>
struct DoesNotHaveAParaboloidReconstruction {
  static constexpr bool value =
      !HasAParaboloidReconstruction<ReconstructionType>::value;
};

template <class ReconstructionType>
struct IsParaboloidReconstruction {
  static constexpr bool value =
      is_paraboloid_reconstruction<ReconstructionType>::value;
};

template <class ReconstructionType>
struct IsNotAParaboloidReconstruction {
  static constexpr bool value =
      !IsParaboloidReconstruction<ReconstructionType>::value;
};

template <class ReconstructionType>
struct IsPlanarSeparatorPathGroup {
  static constexpr bool value =
      has_planar_separator_path_group<ReconstructionType>::value;
};

template <class ReconstructionType>
struct IsNotAPlanarSeparatorPathGroup {
  static constexpr bool value =
      !IsPlanarSeparatorPathGroup<ReconstructionType>::value;
};

template <class ReconstructionType>
struct IsPlanarSeparator {
  static constexpr bool value = is_planar_separator<ReconstructionType>::value;
};

template <class ReconstructionType>
struct IsNotAPlanarSeparator {
  static constexpr bool value = !IsPlanarSeparator<ReconstructionType>::value;
};

template <class ReconstructionType>
struct HasALocalizer {
  static constexpr bool value = has_localizer<ReconstructionType>::value;
};

template <class ReconstructionType>
struct DoesNotHaveALocalizer {
  static constexpr bool value = !HasALocalizer<ReconstructionType>::value;
};

template <class ReconstructionType>
struct HasASeparator {
  static constexpr bool value = has_separator<ReconstructionType>::value;
};

template <class ReconstructionType>
struct DoesNotHaveASeparator {
  static constexpr bool value = !HasASeparator<ReconstructionType>::value;
};

template <class ReconstructionType>
struct HasLocalizer_AND_HasSeparator {
  static constexpr bool value = (HasALocalizer<ReconstructionType>::value &&
                                 HasASeparator<ReconstructionType>::value);
};

template <class ReconstructionType>
struct HasLocalizer_OR_HasSeparator {
  static constexpr bool value = (HasALocalizer<ReconstructionType>::value ||
                                 HasASeparator<ReconstructionType>::value);
};

template <class ReconstructionType>
struct HasLocalizer_XOR_HasSeparator {
  static constexpr bool value =
      (HasLocalizer_OR_HasSeparator<ReconstructionType>::value &&
       !HasLocalizer_AND_HasSeparator<ReconstructionType>::value);
};

template <class ReconstructionType>
struct HasAReconstructionLink {
  static constexpr bool value =
      is_reconstruction_link<ReconstructionType>::value;
};

template <class ReconstructionType>
struct DoesNotHaveAReconstructionLink {
  static constexpr bool value =
      !HasAReconstructionLink<ReconstructionType>::value;
};

template <class ReturnType>
struct HasANestedType {
  static constexpr bool value = is_nested_moments<ReturnType>::value;
};

template <class ReturnType>
struct DoesNotHaveANestedType {
  static constexpr bool value = !HasANestedType<ReturnType>::value;
};

template <class ReturnType>
struct HasACollection {
  static constexpr bool value = is_moments_collection<ReturnType>::value;
};

template <class ReturnType>
struct DoesNotHaveACollection {
  static constexpr bool value = !HasACollection<ReturnType>::value;
};

template <class ReturnType>
struct IsOnlyVolume {
  static constexpr bool value = is_moments_volume<ReturnType>::value;
};

template <class ReturnType>
struct IsNotVolume {
  static constexpr bool value = !IsOnlyVolume<ReturnType>::value;
};
}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_GENERAL_CLASS_CLASSIFICATIONS_H_
