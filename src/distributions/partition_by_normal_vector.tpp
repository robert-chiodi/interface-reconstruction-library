// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DISTRIBUTIONS_PARTITION_BY_NORMAL_VECTOR_TPP_
#define SRC_DISTRIBUTIONS_PARTITION_BY_NORMAL_VECTOR_TPP_

namespace IRL {

template <class MomentsContainerType>
PartitionByNormal<MomentsContainerType>::PartitionByNormal(void)
    : list_ptr_m(nullptr) {}

template <class MomentsContainerType>
PartitionByNormal<MomentsContainerType>::PartitionByNormal(
    const MomentsContainerType* a_polygon_container_ptr)
    : list_ptr_m(a_polygon_container_ptr) {}

template <class MomentsContainerType>
void PartitionByNormal<MomentsContainerType>::setup(
    const MomentsContainerType* a_polygon_container_ptr) {
  list_ptr_m = a_polygon_container_ptr;
  this->setup();
}

template <class MomentsContainerType>
void PartitionByNormal<MomentsContainerType>::setup(void) {
  assert(list_ptr_m != nullptr);
  old_partition_normals_m[0] = 0.0;
  old_partition_normals_m[1] = 0.0;
  // Select initial first normal as one with largest corresponding Volume
  partitioned_objects_m[0] = this->getObjectWithMostVolume();
  Normal normal_for_object_with_most_volume = partitioned_objects_m[0].normal();
  normal_for_object_with_most_volume.normalize();
  partitioned_objects_m[1] =
      this->getObjectMostDifferent(normal_for_object_with_most_volume);
}

template <class MomentsContainerType>
std::array<typename MomentsContainerType::contained_type, 2>
PartitionByNormal<MomentsContainerType>::getPartitionedObjects(void) {
  return partitioned_objects_m;
}

template <class MomentsContainerType>
bool PartitionByNormal<MomentsContainerType>::isDone(void) {
  std::array<Normal, 2> partitioned_objects_normalized_normals;
  partitioned_objects_normalized_normals[0] = partitioned_objects_m[0].normal();
  partitioned_objects_normalized_normals[1] = partitioned_objects_m[1].normal();
  partitioned_objects_normalized_normals[0].normalize();
  partitioned_objects_normalized_normals[1].normalize();
  return partitioned_objects_normalized_normals[0] ==
             old_partition_normals_m[0] &&
         partitioned_objects_normalized_normals[1] ==
             old_partition_normals_m[1];
}

template <class MomentsContainerType>
bool PartitionByNormal<MomentsContainerType>::iterationTooHigh(
    const UnsignedIndex_t a_iteration_number) {
  return a_iteration_number > max_iteration_number;
}

template <class MomentsContainerType>
void PartitionByNormal<MomentsContainerType>::setupNextIteration(void) {
  for (UnsignedIndex_t n = 0; n < 2; ++n) {
    old_partition_normals_m[n] = partitioned_objects_m[n].normal();
    old_partition_normals_m[n].normalize();
    partitioned_objects_m[n] = 0.0;
  }
}

template <class MomentsContainerType>
UnsignedIndex_t PartitionByNormal<MomentsContainerType>::findCorrectPartition(
    const typename MomentsContainerType::contained_type& a_element) {
  Normal element_normal = a_element.normal();
  element_normal.normalize();
  return element_normal * old_partition_normals_m[0] >=
                 element_normal * old_partition_normals_m[1]
             ? 0
             : 1;
}

template <class MomentsContainerType>
void PartitionByNormal<MomentsContainerType>::addElementToPartition(
    const UnsignedIndex_t a_partition,
    const typename MomentsContainerType::contained_type& a_element) {
  partitioned_objects_m[a_partition] += a_element;
}

template <class MomentsContainerType>
typename PartitionByNormal<MomentsContainerType>::const_iterator
PartitionByNormal<MomentsContainerType>::begin(void) noexcept {
  return list_ptr_m->begin();
}
template <class MomentsContainerType>
typename PartitionByNormal<MomentsContainerType>::const_iterator
PartitionByNormal<MomentsContainerType>::begin(void) const noexcept {
  return this->cbegin();
}
template <class MomentsContainerType>
typename PartitionByNormal<MomentsContainerType>::const_iterator
PartitionByNormal<MomentsContainerType>::end(void) const noexcept {
  return this->cend();
}
template <class MomentsContainerType>
typename PartitionByNormal<MomentsContainerType>::const_iterator
PartitionByNormal<MomentsContainerType>::cbegin(void) const noexcept {
  return list_ptr_m->cbegin();
}
template <class MomentsContainerType>
typename PartitionByNormal<MomentsContainerType>::const_iterator
PartitionByNormal<MomentsContainerType>::end(void) noexcept {
  return list_ptr_m->end();
}
template <class MomentsContainerType>
typename PartitionByNormal<MomentsContainerType>::const_iterator
PartitionByNormal<MomentsContainerType>::cend(void) const noexcept {
  return list_ptr_m->cend();
}

template <class MomentsContainerType>
typename MomentsContainerType::contained_type
PartitionByNormal<MomentsContainerType>::getObjectWithMostVolume(void) {
  double max_volume = -DBL_MAX;
  typename MomentsContainerType::contained_type object_with_most_volume;
  for (const auto& element : (*list_ptr_m)) {
    if (element.volumeMoments().volume() > max_volume) {
      max_volume = element.volumeMoments().volume();
      object_with_most_volume = element;
    }
  }
  return object_with_most_volume;
}

template <class MomentsContainerType>
typename MomentsContainerType::contained_type
PartitionByNormal<MomentsContainerType>::getObjectMostDifferent(
    const Normal& a_normal) {
  double max_difference = DBL_MAX;
  typename MomentsContainerType::contained_type object_with_most_difference;
  for (const auto& element : (*list_ptr_m)) {
    Normal normalized_element_normal = element.normal();
    normalized_element_normal.normalize();
    double difference = a_normal * (normalized_element_normal);
    if (difference < max_difference &&
        magnitude(normalized_element_normal) > 0.9999) {
      max_difference = difference;
      object_with_most_difference = element;
    }
  }
  return object_with_most_difference;
}

}  // namespace IRL

#endif  // SRC_DISTRIBUTIONS_PARTITION_BY_NORMAL_VECTOR_TPP_
