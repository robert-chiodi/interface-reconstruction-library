// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_LISTED_VOLUME_MOMENTS_H_
#define SRC_MOMENTS_LISTED_VOLUME_MOMENTS_H_

#include <iomanip>
#include <iostream>

#include "src/data_structures/lister.h"
#include "src/parameters/defined_types.h"

namespace IRL {

/// \brief VolumeMoments wrapper for Lister class
/// to add ability to normalize.
template <class VolumeMomentsType>
class ListedVolumeMoments {
 public:
  using contained_type = VolumeMomentsType;
  using iterator = typename Lister<VolumeMomentsType>::iterator;
  using const_iterator = typename Lister<VolumeMomentsType>::const_iterator;

  /// \brief Default construcor.
  ListedVolumeMoments(void) = default;

  /// \brief Get size of the collection.
  UnsignedIndex_t size(void) const;

  /// \brief This will self-expand to prevent itself from
  /// accessing out of bounds memory.
  contained_type& operator[](const UnsignedIndex_t a_index);

  /// \brief Const version for access to object in collection.
  const contained_type& operator[](const UnsignedIndex_t a_index) const;

  /// \brief Set all moments in the list equal to a_value.
  ListedVolumeMoments& operator=(const double a_value);

  /// \brief Normalize entire vector by volume.
  void normalizeByVolume(void);

  /// \brief Normalize entire vector by volume.
  void multiplyByVolume(void);

  /// \brief Overload operator+= to accumulate moments with the same index, and
  /// extend to match `a_rhs` length.
  ListedVolumeMoments& operator+=(const ListedVolumeMoments& a_rhs);

  /// \brief The operator += will be used to push_back
  /// the object `a_rhs` in the current collection in Lister.
  ListedVolumeMoments& operator+=(const VolumeMomentsType& a_rhs);

  /// \brief Empty the container.
  void clear(void);

  /// \brief Erase an object from the container at the index.
  void erase(const UnsignedIndex_t a_index);

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Default destrucor.
  ~ListedVolumeMoments(void) = default;

 private:
  Lister<VolumeMomentsType> list_m;
};

//******************************************************************* //
//     Inlined function definitions placed below this.
//******************************************************************* //
template <class VolumeMomentsType>
inline std::ostream& operator<<(
    std::ostream& out,
    const ListedVolumeMoments<VolumeMomentsType>& a_listed_volume_moments);

}  // namespace IRL

#include "src/moments/listed_volume_moments.tpp"

#endif  // SRC_MOMENTS_LISTED_VOLUME_MOMENTS_H_
