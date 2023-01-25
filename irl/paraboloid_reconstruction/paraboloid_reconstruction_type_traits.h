// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_PARABOLOID_RECONSTRUCTION_TYPE_TRAITS_H_
#define IRL_PARABOLOID_RECONSTRUCTION_PARABOLOID_RECONSTRUCTION_TYPE_TRAITS_H_

#include "irl/paraboloid_reconstruction/paraboloid.h"
#include "irl/planar_reconstruction/planar_reconstruction_type_traits.h"

namespace IRL {
//******************************************************************* //
//   Is a Paraboloid based reconstruction                       //
//******************************************************************* //
template <class C>
struct is_paraboloid_reconstruction : std::false_type {};

template <class C>
struct is_paraboloid_reconstruction<const C> : is_paraboloid_reconstruction<C> {
};

template <class ScalarType>
struct is_paraboloid_reconstruction<ParaboloidBase<ScalarType>>
    : std::true_type {};

//******************************************************************* //
//   Contains a Paraboloid based reconstruction                       //
//******************************************************************* //
template <class C>
struct has_paraboloid_reconstruction : std::false_type {};

template <class C>
struct has_paraboloid_reconstruction<const C>
    : has_paraboloid_reconstruction<C> {};

template <class ScalarType>
struct has_paraboloid_reconstruction<ParaboloidBase<ScalarType>>
    : std::true_type {};

template <class ScalarType>
struct has_paraboloid_reconstruction<LocalizedParaboloid<ScalarType>>
    : std::true_type {};

template <class ScalarType>
struct has_paraboloid_reconstruction<LocalizedParaboloidLink<ScalarType>>
    : std::true_type {};

//********************************************************************* //
//  Contains a localizer, base from planar_reconstruction_type_traits.h //
//********************************************************************* //
template <class ScalarType>
struct has_localizer<LocalizedParaboloid<ScalarType>> : std::true_type {};

template <class ScalarType>
struct has_localizer<LocalizedParaboloidLink<ScalarType>> : std::true_type {};

//******************************************************************* //
//        Linked Planar Reconstructions
//        True for anything that inherits from UnDirectedGraphNode<Self>
//******************************************************************* //
template <class ScalarType>
struct is_reconstruction_link<LocalizedParaboloidLink<ScalarType>>
    : std::true_type {};

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTION_PARABOLOID_RECONSTRUCTION_TYPE_TRAITS_H_
