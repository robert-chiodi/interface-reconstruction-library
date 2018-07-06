// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_HELPERS_SERIALIZER_TPP_
#define SRC_HELPERS_SERIALIZER_TPP_

namespace IRL {

template <class ObjectType, class ContainerType>
void serializeAndPack(const ObjectType& a_object, ContainerType* a_container) {
  a_object.serialize(a_container);
}

template <class ObjectType, class ContainerType>
void unpackAndStore(ObjectType* a_object, ContainerType* a_container) {
  assert(a_container != nullptr);
  a_object->unpackSerialized(a_container);
}

}  // namespace IRL

#endif  // SRC_HELPERS_SERIALIZER_TPP_
