// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_CELL_COLLECTION_H_
#define SRC_MOMENTS_CELL_COLLECTION_H_

#include "src/data_structures/collection.h"
#include "src/parameters/defined_types.h"

namespace IRL {
/// \brief Class that contains a list of cells and coupled to moments.
template <class CellGroupType>
class CellCollection {
 public:
  using iterator = typename Collection<CellGroupType>::iterator;
  using const_iterator = typename Collection<CellGroupType>::const_iterator;
  using contained_type = CellGroupType;

  /// \brief Default constructor.
  CellCollection(void) = default;

  const typename CellGroupType::cell_type& getCell(
      const UnsignedIndex_t a_index) const;

  const typename CellGroupType::contained_type& getStoredMoments(
      const UnsignedIndex_t a_index) const;

  /// \brief Get size of the collection.
  UnsignedIndex_t size(void) const;

  /// \brief Set number of cells that are going to be involved in the
  /// collection.
  void resize(const UnsignedIndex_t a_size);

  /// \brief Add member to end of collection vector.
  void push_back(const CellGroupType& a_object);

  /// \brief Empty the container.
  void clear(void);

  /// \brief Get a member from the collection.
  CellGroupType& operator[](const UnsignedIndex_t a_index);

  /// \brief Const version of getMember.
  const CellGroupType& operator[](const UnsignedIndex_t a_index) const;

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Default destructor.
  ~CellCollection(void) = default;

 private:
  void checkIndex(const UnsignedIndex_t a_index) const;

  /// \brief Collection of CellGroupedMoments
  Collection<CellGroupType> collection_m;
};

}  // namespace IRL

#include "src/moments/cell_collection.tpp"

#endif  // SRC_MOMENTS_CELL_COLLECTION_H_
