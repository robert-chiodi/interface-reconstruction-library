// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DISTRIBUTIONS_K_MEANS_H_
#define SRC_DISTRIBUTIONS_K_MEANS_H_

#include <cassert>

#include "src/parameters/defined_types.h"

namespace IRL {

/// \brief A class that executes Kmeans when
/// provided an approprimate DrivingClass.
///
/// This function is used to drive Kmeans in order
/// to partition a set into multiple sets. This is
/// done in a very broad way which allows it to perform
/// this partitioning and the effect of the partitioning
/// in many ways, dictated by the driving class.
///
/// Requirements for DrivingClass:
/// - `bool isDone(void)` : A method that determines if
/// the KMeans routine should be stopped and the initial set
/// has been sufficiently partitioned.
/// -`void setupNextIteration(void)` : A method to prepare
/// everything necessary in the DrivingClass object for the
/// next iteration/partitioning.
/// - `const_iterator begin()` : Method that returns a const
/// iterator to the start of the initial set being partitioned.
/// - `const_iterator end()` : Method that returns a const
/// iterator to the end of the initial set being partitioned.
/// -`int findCorrectPartition(const ElementType&)` : A method that
/// takes an element from the initial set in the DrivingClass
/// object and returns the partition in DrivingClass that it should
/// belong to.
/// -`void addElementToPartition(const int, const ElementType&)` :
/// Given the partition the element should belong to (determined by the
/// `findCorrectPartition` method, the element is added to the partition.
/// Can be used to directly update the consequence instead, such as
/// updating a sum instead of directly constructing sets representing
/// the partition.
/// -`bool iterationTooHigh(const UnsignedIndex_t)` : A method that takes
/// the number of iterations and returns a bool
/// whether the maximum number of allowable iterations
/// has been exceeded.

struct KMeans {
  template <class DrivingClass>
  /// \brief Partition into KMeans according to DrivingClass.
  static UnsignedIndex_t partition(DrivingClass* a_driving_ptr);
};

}  // namespace IRL

#include "src/distributions/k_means.tpp"

#endif  // SRC_DISTRIBUTIONS_K_MEANS_H_
