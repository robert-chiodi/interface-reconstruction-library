// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DISTRIBUTIONS_PARTITION_BY_NORMAL_VECTOR_H_
#define SRC_DISTRIBUTIONS_PARTITION_BY_NORMAL_VECTOR_H_

#include <cstring>
#include <utility>

#include "src/distributions/k_means.h"
#include "src/geometry/general/normal.h"
#include "src/parameters/defined_types.h"

namespace IRL {

/// \brief This is a class that takes a list of objects
/// and separates the list into two partitions, which
/// are represented by the sums of the objects.
///
/// This class is used to split a list of objects
/// into two separate groups depending on normal vectors.
/// The normal vectors contained as part of the objects
/// in the list are assumed to not necessarily be normalized,
/// and normalization is done before they are used. Instead of saving the
/// elements in each partition, the sum of the objects (through operator +=) is
/// computed for each partition. This is done by providing this class (once
/// setup) to the KMeans class, where KMeans::solve() is then used.
///
///
///
/// Requirements for ContainerType:
/// - Const_Iterator that allows traversing the entire initial set
/// held in list_ptr_m.
/// - `const Normal& normal(void)` : A method that returns a Normal
/// (or const reference to one) to be used for comparison. This
/// Normal does not need to be normalized, will be normalized during
/// partitioning.
/// `Object& overload+=(const Object&)` : A method that overloads +=
/// to enable summing of the elements into their partitions.
template <class MomentsContainerType>
class PartitionByNormal {
  friend KMeans;

  using iterator = typename MomentsContainerType::iterator;
  using const_iterator = typename MomentsContainerType::const_iterator;
  static constexpr UnsignedIndex_t max_iteration_number = 40;

 public:
  /// \brief Default constructor.
  PartitionByNormal(void);

  /// \brief Constructor that takes pointer to list of polygons to be
  /// partitioned.
  explicit PartitionByNormal(
      const MomentsContainerType* a_polygon_container_ptr);

  /// \brief Set pointer, then setup the system before passing to KMeans.
  void setup(const MomentsContainerType* a_polygon_container_ptr);

  /// \brief Setup the system before passing to KMeans.
  void setup(void);

  /// \brief Write partitioned objects to a provided buffer.
  std::array<typename MomentsContainerType::contained_type, 2>
  getPartitionedObjects(void);

  /// \brief Default constructor.
  ~PartitionByNormal(void) = default;

 private:
  /// \brief Check if partitioning isdone by seeing if the found normals
  /// did not change during iteration.
  bool isDone(void);

  /// \brief Return bool indicating if max iterations has
  /// been surpassed.
  bool iterationTooHigh(const UnsignedIndex_t a_iteration_number);

  /// \brief Setup the next iteration by updating old normals
  /// to previous iteration value and zeroing partitioned moments.
  void setupNextIteration(void);

  /// \brief Find which partition the normal belongs in.
  UnsignedIndex_t findCorrectPartition(
      const typename MomentsContainerType::contained_type& a_element);

  /// \brief Add the element to the indicated partition.
  void addElementToPartition(
      const UnsignedIndex_t a_partition,
      const typename MomentsContainerType::contained_type& a_element);

  const_iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  const_iterator end(void) noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Find and return the object in the list with the most associated
  /// volume.
  typename MomentsContainerType::contained_type getObjectWithMostVolume(void);

  /// \brief Find and return the object in the list with the normal most
  /// different from a_normal.
  typename MomentsContainerType::contained_type getObjectMostDifferent(
      const Normal& a_normal);

  const MomentsContainerType* list_ptr_m;
  std::array<typename MomentsContainerType::contained_type, 2>
      partitioned_objects_m;
  std::array<Normal, 2> old_partition_normals_m;
};

}  // namespace IRL

#include "src/distributions/partition_by_normal_vector.tpp"

#endif  // SRC_DISTRIBUTIONS_PARTITION_BY_NORMAL_VECTOR_H_
