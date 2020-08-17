// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_LISTER_TPP_
#define SRC_DATA_STRUCTURES_LISTER_TPP_

namespace IRL {

template <class ObjectType>
Lister<ObjectType>& Lister<ObjectType>::operator+=(const Lister& a_rhs) {
  const UnsignedIndex_t old_size = this->size();
  this->resize(old_size + a_rhs.size());
  for (UnsignedIndex_t i = old_size; i < (old_size + a_rhs.size()); ++i) {
    (*this)[i] = a_rhs[i - old_size];
  }
  return (*this);
}

template <class ObjectType>
Lister<ObjectType>& Lister<ObjectType>::operator+=(const ObjectType& a_rhs) {
  this->push_back(a_rhs);
  return (*this);
}

}  // namespace IRL

#endif  // SRC_DATA_STRUCTURES_LISTER_TPP_
