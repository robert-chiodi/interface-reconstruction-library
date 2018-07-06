// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_HELPERS_SERIALIZER_H_
#define SRC_HELPERS_SERIALIZER_H_

#include <cassert>

namespace IRL {

// Bit conversion is written specifically for homogeneous architectures!! Little
// endian -> Little endian or Big endian -> Big endian

template <class ObjectType, class ContainerType>
void serializeAndPack(const ObjectType& a_object, ContainerType* a_container);

template <class ObjectType, class ContainerType>
void unpackAndStore(ObjectType* a_object, ContainerType* a_container);

}  // namespace IRL

#include "src/helpers/serializer.tpp"

#endif  // SRC_HELPERS_SERIALIZER_H_
