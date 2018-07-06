// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_TAGGED_ACCUMULATED_VOLUME_MOMENTS_TPP_
#define SRC_MOMENTS_TAGGED_ACCUMULATED_VOLUME_MOMENTS_TPP_

#include <cassert>
#include <iostream>

namespace IRL {

template <class VolumeMomentsType>
TaggedAccumulatedVolumeMoments<
    VolumeMomentsType>::TaggedAccumulatedVolumeMoments(void)
    : accumulated_moments_m(), tag_to_vector_index_m() {
  // tag_to_vector_index_m.reserve(initial_bucket_count);
}

template <class VolumeMomentsType>
TaggedAccumulatedVolumeMoments<VolumeMomentsType>::
    TaggedAccumulatedVolumeMoments(
        const TaggedAccumulatedVolumeMoments &a_other) noexcept
    : accumulated_moments_m(a_other.accumulated_moments_m),
      tag_to_vector_index_m(a_other.tag_to_vector_index_m) {}

template <class VolumeMomentsType>
TaggedAccumulatedVolumeMoments<VolumeMomentsType>
    &TaggedAccumulatedVolumeMoments<VolumeMomentsType>::operator=(
        const TaggedAccumulatedVolumeMoments &a_other) noexcept {
  if (this != &a_other) {
    accumulated_moments_m = a_other.accumulated_moments_m;
    tag_to_vector_index_m = a_other.tag_to_vector_index_m;
  }
  return *this;
}

template <class VolumeMomentsType>
VolumeMomentsType &TaggedAccumulatedVolumeMoments<VolumeMomentsType>::
operator[](const UnsignedIndex_t a_tag) {
  const auto insert_pair = tag_to_vector_index_m.insert({a_tag, this->size()});
  if (insert_pair.second) {  // Successful insertion
    accumulated_moments_m.emplace_back(
        TaggedMoments(VolumeMomentsType(), a_tag));
    return accumulated_moments_m[this->size() - 1].volume_moments_m;
  }
  return accumulated_moments_m[(*insert_pair.first).second].volume_moments_m;
}

template <class VolumeMomentsType>
const VolumeMomentsType &TaggedAccumulatedVolumeMoments<VolumeMomentsType>::
operator[](const UnsignedIndex_t a_tag) const {
  assert(this->isTagKnown(a_tag));
  return accumulated_moments_m[tag_to_vector_index_m.at(a_tag)]
      .volume_moments_m;
}

template <class VolumeMomentsType>
VolumeMomentsType &
TaggedAccumulatedVolumeMoments<VolumeMomentsType>::getMomentsForIndex(
    const UnsignedIndex_t a_index) {
  assert(a_index < this->size());
  return accumulated_moments_m[a_index].volume_moments_m;
}

template <class VolumeMomentsType>
const VolumeMomentsType &
TaggedAccumulatedVolumeMoments<VolumeMomentsType>::getMomentsForIndex(
    const UnsignedIndex_t a_index) const {
  assert(a_index < this->size());
  return accumulated_moments_m[a_index].volume_moments_m;
}

template <class VolumeMomentsType>
UnsignedIndex_t
TaggedAccumulatedVolumeMoments<VolumeMomentsType>::getTagForIndex(
    const UnsignedIndex_t a_index) const {
  assert(a_index < this->size());
  return accumulated_moments_m[a_index].tag_m;
}

template <class VolumeMomentsType>
TaggedAccumulatedVolumeMoments<VolumeMomentsType>
    &TaggedAccumulatedVolumeMoments<VolumeMomentsType>::operator+=(
        const TaggedAccumulatedVolumeMoments &a_other) {
  for (UnsignedIndex_t elem = 0; elem < a_other.size(); ++elem) {
    UnsignedIndex_t tag_for_others_element = a_other.getTagForIndex(elem);
    // Try to insert, if fails, accumulate. If successful, push back onto
    // vector.
    const auto insert_pair =
        tag_to_vector_index_m.insert({tag_for_others_element, this->size()});
    if (insert_pair.second) {  // Successful insertion
      accumulated_moments_m.emplace_back(TaggedMoments(
          a_other.getMomentsForIndex(elem), tag_for_others_element));
    } else {  // Tag already existed
      accumulated_moments_m[(*insert_pair.first).second].volume_moments_m +=
          a_other.getMomentsForIndex(elem);
    }
  }
  return (*this);
}

template <class VolumeMomentsType>
bool TaggedAccumulatedVolumeMoments<VolumeMomentsType>::isTagKnown(
    const UnsignedIndex_t a_tag) const {
  return tag_to_vector_index_m.find(a_tag) != tag_to_vector_index_m.end();
}

template <class VolumeMomentsType>
bool TaggedAccumulatedVolumeMoments<VolumeMomentsType>::isTagNew(
    const UnsignedIndex_t a_tag) const {
  return !this->isTagKnown(a_tag);
}

template <class VolumeMomentsType>
void TaggedAccumulatedVolumeMoments<VolumeMomentsType>::normalizeByVolume(
    void) {
  for (auto &volume_moment : (*this)) {
    volume_moment.volume_moments_m.normalizeByVolume();
  }
}

template <class VolumeMomentsType>
void TaggedAccumulatedVolumeMoments<VolumeMomentsType>::multiplyByVolume(void) {
  for (auto &volume_moment : (*this)) {
    volume_moment.volume_moments_m.multiplyByVolume();
  }
}

template <class VolumeMomentsType>
TaggedAccumulatedVolumeMoments<VolumeMomentsType>
    &TaggedAccumulatedVolumeMoments<VolumeMomentsType>::operator=(
        const double a_value) {
  for (auto &moments : (*this)) {
    moments.volume_moments_m = a_value;
  }
  return (*this);
}

template <class VolumeMomentsType>
UnsignedIndex_t TaggedAccumulatedVolumeMoments<VolumeMomentsType>::size(
    void) const {
  return static_cast<UnsignedIndex_t>(accumulated_moments_m.size());
}

template <class VolumeMomentsType>
void TaggedAccumulatedVolumeMoments<VolumeMomentsType>::addNewTaggedAddress(
    const UnsignedIndex_t a_tag) {
  assert(this->isTagNew(a_tag));
  tag_to_vector_index_m[a_tag] = this->size();
  accumulated_moments_m.resize(this->size() + 1);
  accumulated_moments_m.back().tag_m = a_tag;
}

template <class VolumeMomentsType>
void TaggedAccumulatedVolumeMoments<VolumeMomentsType>::clear(void) {
  accumulated_moments_m.clear();
  tag_to_vector_index_m.clear();
}

template <class VolumeMomentsType>
typename TaggedAccumulatedVolumeMoments<VolumeMomentsType>::iterator
TaggedAccumulatedVolumeMoments<VolumeMomentsType>::begin(void) noexcept {
  return accumulated_moments_m.begin();
}
template <class VolumeMomentsType>
typename TaggedAccumulatedVolumeMoments<VolumeMomentsType>::const_iterator
TaggedAccumulatedVolumeMoments<VolumeMomentsType>::begin(void) const noexcept {
  return this->cbegin();
}
template <class VolumeMomentsType>
typename TaggedAccumulatedVolumeMoments<VolumeMomentsType>::const_iterator
TaggedAccumulatedVolumeMoments<VolumeMomentsType>::cbegin(void) const noexcept {
  return accumulated_moments_m.cbegin();
}
template <class VolumeMomentsType>
typename TaggedAccumulatedVolumeMoments<VolumeMomentsType>::iterator
TaggedAccumulatedVolumeMoments<VolumeMomentsType>::end(void) noexcept {
  return accumulated_moments_m.end();
}
template <class VolumeMomentsType>
typename TaggedAccumulatedVolumeMoments<VolumeMomentsType>::const_iterator
TaggedAccumulatedVolumeMoments<VolumeMomentsType>::end(void) const noexcept {
  return this->cend();
}
template <class VolumeMomentsType>
typename TaggedAccumulatedVolumeMoments<VolumeMomentsType>::const_iterator
TaggedAccumulatedVolumeMoments<VolumeMomentsType>::cend(void) const noexcept {
  return accumulated_moments_m.cend();
}

template <class VolumeMomentsType>
inline std::ostream &operator<<(
    std::ostream &out,
    const TaggedAccumulatedVolumeMoments<VolumeMomentsType> &a_list) {
  out << "There are " << a_list.size() << " elements in the list." << std::endl;
  for (UnsignedIndex_t element = 0; element < a_list.size(); ++element) {
    out << "Tag: " << a_list.getTagForIndex(element)
        << "  Moment: " << a_list.getMomentsForIndex(element) << std::endl;
  }
  return out;
}

template <class VolumeMomentsType>
inline TaggedAccumulatedVolumeMoments<VolumeMomentsType> operator*(
    const TaggedAccumulatedVolumeMoments<VolumeMomentsType> &a_list,
    const double a_multiplier) {
  auto list_to_return = a_list;
  for (UnsignedIndex_t elem = 0; elem < list_to_return.size(); ++elem) {
    list_to_return.getMomentsForIndex(elem) *= a_multiplier;
  }
  return list_to_return;
}

template <class VolumeMomentsType>
inline TaggedAccumulatedVolumeMoments<VolumeMomentsType> operator*(
    const double a_multiplier,
    const TaggedAccumulatedVolumeMoments<VolumeMomentsType> &a_list) {
  return a_list * a_multiplier;
}

}  // namespace IRL

#endif  // SRC_MOMENTS_TAGGED_ACCUMULATED_VOLUME_MOMENTS_TPP_
