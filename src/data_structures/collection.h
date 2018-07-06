// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_COLLECTION_H_
#define SRC_DATA_STRUCTURES_COLLECTION_H_

#include <cassert>
#include <vector>

#include "src/parameters/defined_types.h"

namespace IRL {
/// \brief Essentially a std::vector with
/// assertions on access through operator[] in
/// debug mode.
template <class ObjectType>
class Collection {
 public:
  using iterator = typename std::vector<ObjectType>::iterator;
  using const_iterator = typename std::vector<ObjectType>::const_iterator;
  using contained_type = ObjectType;

  /// \brief Default constructor.
  Collection(void) = default;

  /// \brief Get size of the collection.
  UnsignedIndex_t size(void) const;

  /// \brief Set number of cells that are going to be involved in the
  /// collection.
  void resize(const UnsignedIndex_t a_size);

  /// \brief Add member to end of collection vector.
  void push_back(const ObjectType& a_object);

  /// \brief Construct member at end of collection vector.
  void emplace_back(const ObjectType& a_object);

  /// \brief Empty the container.
  void clear(void);

  /// \brief Erase an object from the container at the index.
  void erase(const UnsignedIndex_t a_index);

  /// \brief Get a member from the collection.
  ObjectType& operator[](const UnsignedIndex_t a_index);

  /// \brief Const version of getMember.
  const ObjectType& operator[](const UnsignedIndex_t a_index) const;

  iterator begin(void) noexcept;
  const_iterator begin(void) const noexcept;
  const_iterator cbegin(void) const noexcept;
  iterator end(void) noexcept;
  const_iterator end(void) const noexcept;
  const_iterator cend(void) const noexcept;

  /// \brief Default destructor.
  ~Collection(void) = default;

 protected:
  void checkIndex(const UnsignedIndex_t a_index) const;

 private:
  std::vector<ObjectType> collection_m;
};
}  // namespace IRL

#include "src/data_structures/collection.tpp"

#endif  // SRC_DATA_STRUCTURES_COLLECTION_H_
