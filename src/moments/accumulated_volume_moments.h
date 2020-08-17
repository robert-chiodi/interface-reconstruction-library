// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_ACCUMULATED_VOLUME_MOMENTS_H_
#define SRC_MOMENTS_ACCUMULATED_VOLUME_MOMENTS_H_

#include "src/data_structures/accumulator.h"
#include "src/parameters/defined_types.h"

namespace IRL {

/// \brief Self-expanding vector of VolumeMoments.
template <class VolumeMomentsType>
class AccumulatedVolumeMoments {
 public:
  using iterator = typename Accumulator<VolumeMomentsType>::iterator;
  using const_iterator =
      typename Accumulator<VolumeMomentsType>::const_iterator;
  using contained_type = VolumeMomentsType;

  /// \brief Default constructor
  AccumulatedVolumeMoments(void) = default;

  /// \brief Get size of the collection.
  UnsignedIndex_t size(void) const;

  /// \brief This will self-expand to prevent itself from
  /// accessing out of bounds memory.
  contained_type& operator[](const UnsignedIndex_t a_index);

  /// \brief Const version for access to object in collection.
  const contained_type& operator[](const UnsignedIndex_t a_index) const;

  /// \brief Set all moments in the list equal to a_value.
  AccumulatedVolumeMoments& operator=(const double a_value);

  /// \brief Normalize entire vector by volume.
  void normalizeByVolume(void);

  /// \brief Normalize entire vector by volume.
  void multiplyByVolume(void);

  /// \brief Overload operator+= to accumulate moments with the same index, and
  /// extend to match `a_rhs` length.
  AccumulatedVolumeMoments& operator+=(const AccumulatedVolumeMoments& a_rhs);

  /// \brief Empty the container.
  void clear(void);

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Default destructor
  ~AccumulatedVolumeMoments(void) = default;

 private:
  Accumulator<VolumeMomentsType> accumulated_moments_m;
};

/// \brief Overload * operator to multiply the contained VolumeMomentsType.
template <class VolumeMomentsType>
inline AccumulatedVolumeMoments<VolumeMomentsType> operator*(
    const AccumulatedVolumeMoments<VolumeMomentsType>& a_list,
    const double a_multiplier);
/// \brief Overload * operator to multiply the contained VolumeMomentsType.
template <class VolumeMomentsType>
inline AccumulatedVolumeMoments<VolumeMomentsType> operator*(
    const double a_multiplier,
    const AccumulatedVolumeMoments<VolumeMomentsType>& a_list);
}  // namespace IRL

#include "src/moments/accumulated_volume_moments.tpp"

#endif  // SRC_MOMENTS_ACCUMULATED_VOLUME_MOMENTS_H_
