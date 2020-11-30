// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_TAGGED_ACCUMULATED_VOLUME_MOMENTS_H_
#define IRL_MOMENTS_TAGGED_ACCUMULATED_VOLUME_MOMENTS_H_

#include <functional>
#include <unordered_map>
#include <utility>

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/unordered_map.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class VolumeMomentsType>
class TaggedAccumulatedVolumeMoments {
  class TaggedMoments {
   public:
    TaggedMoments(void) = default;
    TaggedMoments(const VolumeMomentsType &a_vm, const UnsignedIndex_t a_tag)
        : volume_moments_m{a_vm}, tag_m{a_tag} {}

    VolumeMomentsType volume_moments_m;
    UnsignedIndex_t tag_m;
  };

  static constexpr std::size_t initial_bucket_count = 10;
  using VectorType = SmallVector<TaggedMoments, initial_bucket_count>;
  using iterator = typename VectorType::iterator;
  using const_iterator = typename VectorType::const_iterator;

 public:
  using contained_type = VolumeMomentsType;
  TaggedAccumulatedVolumeMoments(void);

  TaggedAccumulatedVolumeMoments(
      const TaggedAccumulatedVolumeMoments &a_other) noexcept;
  TaggedAccumulatedVolumeMoments &operator=(
      const TaggedAccumulatedVolumeMoments &a_other) noexcept;

  VolumeMomentsType &operator[](const UnsignedIndex_t a_tag);

  const VolumeMomentsType &operator[](const UnsignedIndex_t a_tag) const;

  VolumeMomentsType &getMomentsForIndex(const UnsignedIndex_t a_index);

  const VolumeMomentsType &getMomentsForIndex(
      const UnsignedIndex_t a_index) const;

  UnsignedIndex_t getTagForIndex(const UnsignedIndex_t a_index) const;

  TaggedAccumulatedVolumeMoments &operator+=(
      const TaggedAccumulatedVolumeMoments &a_other);

  bool isTagKnown(const UnsignedIndex_t a_tag) const;

  bool isTagNew(const UnsignedIndex_t a_tag) const;

  /// \brief Normalize entire vector by volume.
  void normalizeByVolume(void);

  /// \brief Normalize entire vector by volume.
  void multiplyByVolume(void);

  /// \brief Set all moments in the list equal to a_value.
  TaggedAccumulatedVolumeMoments &operator=(const double a_value);

  UnsignedIndex_t size(void) const;

  void addNewTaggedAddress(const UnsignedIndex_t a_tag);

  void clear(void);

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  ~TaggedAccumulatedVolumeMoments(void) = default;

 private:
  VectorType accumulated_moments_m;
  IRL::unordered_map<UnsignedIndex_t, UnsignedIndex_t> tag_to_vector_index_m;
};

template <class VolumeMomentsType>
inline std::ostream &operator<<(
    std::ostream &out,
    const TaggedAccumulatedVolumeMoments<VolumeMomentsType> &a_list);

/// \brief Overload * operator to multiply the contained VolumeMomentsType.
template <class VolumeMomentsType>
inline TaggedAccumulatedVolumeMoments<VolumeMomentsType> operator*(
    const TaggedAccumulatedVolumeMoments<VolumeMomentsType> &a_list,
    const double a_multiplier);
/// \brief Overload * operator to multiply the contained VolumeMomentsType.
template <class VolumeMomentsType>
inline TaggedAccumulatedVolumeMoments<VolumeMomentsType> operator*(
    const double a_multiplier,
    const TaggedAccumulatedVolumeMoments<VolumeMomentsType> &a_list);

}  // namespace IRL

#include "irl/moments/tagged_accumulated_volume_moments.tpp"

#endif // IRL_MOMENTS_TAGGED_ACCUMULATED_VOLUME_MOMENTS_H_
