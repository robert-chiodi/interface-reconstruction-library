// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_MOMENTS_TYPE_TRAITS_H_
#define SRC_MOMENTS_MOMENTS_TYPE_TRAITS_H_

#include "src/moments/accumulate_wrapper.h"
#include "src/moments/accumulated_listed_volume_moments.h"
#include "src/moments/accumulated_volume_moments.h"
#include "src/moments/listed_volume_moments.h"
#include "src/moments/separated_volume_moments.h"
#include "src/moments/tagged_accumulated_listed_volume_moments.h"
#include "src/moments/tagged_accumulated_volume_moments.h"
#include "src/moments/volume.h"
#include "src/moments/volume_moments.h"
#include "src/moments/volume_moments_and_normal.h"

namespace IRL {

/// \brief Type trait indicating the moments type is Volume
template <class C>
struct is_moments_volume : std::false_type {};

template <class C>
struct is_moments_volume<const C> : is_moments_volume<C> {};

template <>
struct is_moments_volume<Volume> : std::true_type {};

/// \brief Type trait to allow static checking that an object is a
/// AccumulatedVolumeMoments.
template <class C>
struct is_moments_collection : std::false_type {};

template <class C>
struct is_moments_collection<const C> : is_moments_collection<C> {};

/// \brief Any instantiation of AccumulatedVolumeMoments is a type of
/// volumeMomentsList.
template <class VolumeMomentsType>
struct is_moments_collection<AccumulatedVolumeMoments<VolumeMomentsType>>
    : std::true_type {};

template <class VolumeMomentsType>
struct is_moments_collection<TaggedAccumulatedVolumeMoments<VolumeMomentsType>>
    : std::true_type {};

/// \brief Any instantiation of ListedVolumeMoments is a type of
/// volumeMomentsList.
template <class VolumeMomentsType>
struct is_moments_collection<ListedVolumeMoments<VolumeMomentsType>>
    : std::true_type {};

/// \brief Any instantiation of AccumulatedListedVolumeMoments is a type of
/// volumeMomentsList.
template <class VolumeMomentsType>
struct is_moments_collection<AccumulatedListedVolumeMoments<VolumeMomentsType>>
    : std::true_type {};

/// \brief Any instantiation of AccumulatedListedVolumeMoments is a type of
/// volumeMomentsList.
template <class VolumeMomentsType>
struct is_moments_collection<
    TaggedAccumulatedListedVolumeMoments<VolumeMomentsType>> : std::true_type {
};

/// \brief Type trait to allow static checking that an object is a
/// AccumulatedVolumeMoments.
template <class C>
struct is_nested_moments : is_moments_collection<C> {};

template <class C>
struct is_nested_moments<const C> : is_nested_moments<C> {};

/// \brief Any instantiation of AccumulatedVolumeMoments is a type of
/// volumeMomentsList.
template <class NestedType>
struct is_nested_moments<AccumulateWrapper<NestedType>> : std::true_type {};

/// Marking of if volume
template <class C>
struct is_separated_moments : std::false_type {};

template <class C>
struct is_separated_moments<const C> : is_separated_moments<C> {};

template <class MomentsType>
struct is_separated_moments<SeparatedMoments<MomentsType>> : std::true_type {};

}  // namespace IRL

#endif  // SRC_MOMENTS_MOMENTS_TYPE_TRAITS_H_
