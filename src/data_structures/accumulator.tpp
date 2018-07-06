// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_ACCUMULATOR_TPP_
#define SRC_DATA_STRUCTURES_ACCUMULATOR_TPP_

#include "src/parameters/defined_types.h"

namespace IRL {
template <class ObjectType>
Accumulator<ObjectType>& Accumulator<ObjectType>::operator+=(
    const Accumulator& a_rhs) {
  // This will self expand (*this) to be >= a_rhs.size()
  // Object should be set to default construct with 0's for
  // accumulation to work correctly.
  for (UnsignedIndex_t i = 0; i < a_rhs.size(); ++i) {
    (*this)[i] += a_rhs[i];
  }
  return (*this);
}
}  // namespace IRL

#endif  // SRC_DATA_STRUCTURES_ACCUMULATOR_TPP_
