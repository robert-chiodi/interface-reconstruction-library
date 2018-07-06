// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PLANAR_RECONSTRUCTION_PLANAR_RECONSTRUCTION_TYPE_TRAITS_H_
#define SRC_PLANAR_RECONSTRUCTION_PLANAR_RECONSTRUCTION_TYPE_TRAITS_H_

#include <src/planar_reconstruction/restricted_masked_localized_separator_link.h>
#include "src/graphs/un_directed_graph_node.h"
#include "src/planar_reconstruction/localized_separator.h"
#include "src/planar_reconstruction/localized_separator_group_link.h"
#include "src/planar_reconstruction/localized_separator_link.h"
#include "src/planar_reconstruction/localizer_link.h"
#include "src/planar_reconstruction/localizer_link_from_localized_separator_link.h"
#include "src/planar_reconstruction/planar_localizer.h"
#include "src/planar_reconstruction/masked_localized_separator_link.h"
#include "src/planar_reconstruction/null_reconstruction.h"
#include "src/planar_reconstruction/planar_separator_path.h"
#include "src/planar_reconstruction/planar_separator_path_group.h"
#include "src/planar_reconstruction/restricted_localizer_link_from_localized_separator_link.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {
//******************************************************************* //
//   Is a NullReconstruction                                          //
//******************************************************************* //
template <class C>
struct has_null_reconstruction : std::false_type {};

template <class C>
struct has_null_reconstruction<const C> : has_null_reconstruction<C> {};

template <>
struct has_null_reconstruction<NullReconstruction> : std::true_type {};

//******************************************************************* //
//   Is a PlanarSeparatorPathGroup                                    //
//******************************************************************* //
template <class C>
struct has_planar_separator_path_group : std::false_type {};

template <class C>
struct has_planar_separator_path_group<const C> : has_planar_separator_path_group<C> {};

template <>
struct has_planar_separator_path_group<PlanarSeparatorPathGroup> : std::true_type {};

//******************************************************************* //
//    Is a PlanarSeparator                                            //
//******************************************************************* //
template <class C>
struct is_planar_separator : std::false_type {};

template <class C>
struct is_planar_separator<const C> : is_planar_separator<C> {};

template <>
struct is_planar_separator<PlanarSeparator> : std::true_type {};

//******************************************************************* //
//        Planar Reconstructions that contain a Localizer
//******************************************************************* //
template <class C>
struct has_localizer : std::false_type {};

template <class C>
struct has_localizer<const C> : has_localizer<C> {};

template <>
struct has_localizer<PlanarLocalizer> : std::true_type {};

template <>
struct has_localizer<LocalizedSeparator> : std::true_type {};

template <>
struct has_localizer<LocalizerLink> : std::true_type {};

template <>
struct has_localizer<LocalizedSeparatorLink> : std::true_type {};

template <>
struct has_localizer<MaskedLocalizedSeparatorLink> : std::true_type {};

template <>
struct has_localizer<RestrictedMaskedLocalizedSeparatorLink> : std::true_type {
};

template <>
struct has_localizer<LocalizerLinkFromLocalizedSeparatorLink> : std::true_type {
};

template <>
struct has_localizer<RestrictedLocalizerLinkFromLocalizedSeparatorLink>
    : std::true_type {};

template <>
struct has_localizer<LocalizedSeparatorGroup> : std::true_type {};

template <>
struct has_localizer<LocalizedSeparatorGroupLink> : std::true_type {};

//******************************************************************* //
//        Planar Reconstructions that contain a Separator
//******************************************************************* //
template <class C>
struct has_separator : std::false_type {};

template <class C>
struct has_separator<const C> : has_separator<C> {};

template <>
struct has_separator<PlanarSeparator> : std::true_type {};

template <>
struct has_separator<LocalizedSeparator> : std::true_type {};

template <>
struct has_separator<LocalizedSeparatorLink> : std::true_type {};

template <>
struct has_separator<MaskedLocalizedSeparatorLink> : std::true_type {};

template <>
struct has_separator<RestrictedMaskedLocalizedSeparatorLink> : std::true_type {
};

template <>
struct has_separator<PlanarSeparatorPath> : std::true_type {};

template <>
struct has_separator<LocalizedSeparatorGroup> : std::true_type {};

template <>
struct has_separator<LocalizedSeparatorGroupLink> : std::true_type {};

//******************************************************************* //
//        Linked Planar Reconstructions
//        True for anything that inherits from UnDirectedGraphNode<Self>
//******************************************************************* //
template <class C>
struct is_reconstruction_link : std::false_type {};

template <class C>
struct is_reconstruction_link<const C> : is_reconstruction_link<C> {};

template <>
struct is_reconstruction_link<LocalizerLink> : std::true_type {};

template <>
struct is_reconstruction_link<LocalizedSeparatorLink> : std::true_type {};

template <>
struct is_reconstruction_link<MaskedLocalizedSeparatorLink> : std::true_type {};

template <>
struct is_reconstruction_link<RestrictedMaskedLocalizedSeparatorLink>
    : std::true_type {};

template <>
struct is_reconstruction_link<LocalizerLinkFromLocalizedSeparatorLink>
    : std::true_type {};

template <>
struct is_reconstruction_link<RestrictedLocalizerLinkFromLocalizedSeparatorLink>
    : std::true_type {};

template <>
struct is_reconstruction_link<PlanarSeparatorPath> : std::true_type {};

template <>
struct is_reconstruction_link<LocalizedSeparatorGroupLink> : std::true_type {};

}  // namespace IRL

#endif  // SRC_PLANAR_RECONSTRUCTION_PLANAR_RECONSTRUCTION_TYPE_TRAITS_H_
