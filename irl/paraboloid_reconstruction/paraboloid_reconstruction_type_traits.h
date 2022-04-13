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

#include "irl/paraboloid_reconstruction/aligned_paraboloid.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

namespace IRL {
//******************************************************************* //
//   Is a Paraboloid based reconstruction                             //
//******************************************************************* //
template <class C>
struct has_paraboloid_reconstruction : std::false_type {};

template <class C>
struct has_paraboloid_reconstruction<const C>
    : has_paraboloid_reconstruction<C> {};

template <>
struct has_paraboloid_reconstruction<AlignedParaboloid> : std::true_type {};

template <>
struct has_paraboloid_reconstruction<Paraboloid> : std::true_type {};

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTION_PARABOLOID_RECONSTRUCTION_TYPE_TRAITS_H_
