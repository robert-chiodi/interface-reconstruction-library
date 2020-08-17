// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DISTRIBUTIONS_K_MEANS_TPP_
#define SRC_DISTRIBUTIONS_K_MEANS_TPP_

namespace IRL {

template <class DrivingClass>
/// \brief Partition into KMeans according to DrivingClass.

UnsignedIndex_t KMeans::partition(DrivingClass* a_driving_ptr) {
  assert(a_driving_ptr != nullptr);
  UnsignedIndex_t iteration = 0;
  while (!a_driving_ptr->isDone()) {
    a_driving_ptr->setupNextIteration();
    for (const auto& element : (*a_driving_ptr)) {
      auto belonging_partition = a_driving_ptr->findCorrectPartition(element);
      a_driving_ptr->addElementToPartition(belonging_partition, element);
    }
    ++iteration;
    if (a_driving_ptr->iterationTooHigh(iteration)) {
      break;
    }
  }
  return iteration;
}

}  // namespace IRL

#endif  // SRC_DISTRIBUTIONS_K_MEANS_TPP_
