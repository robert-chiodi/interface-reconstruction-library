// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_SELF_EXPANDING_COLLECTION_TPP_
#define SRC_DATA_STRUCTURES_SELF_EXPANDING_COLLECTION_TPP_

namespace IRL {

template <class ObjectType>
ObjectType& SelfExpandingCollection<ObjectType>::operator[](
    const UnsignedIndex_t a_index) {
  if (a_index + 1 > this->size()) {
    this->resize(a_index + 1);
  }
  return Collection<ObjectType>::operator[](a_index);
}

template <class ObjectType>
const ObjectType& SelfExpandingCollection<ObjectType>::operator[](
    const UnsignedIndex_t a_index) const {
  this->checkIndex(a_index);
  return Collection<ObjectType>::operator[](a_index);
}

}  // namespace IRL

#endif  // SRC_DATA_STRUCTURES_SELF_EXPANDING_COLLECTION_TPP_
