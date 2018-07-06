// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_SELF_EXPANDING_COLLECTION_H_
#define SRC_DATA_STRUCTURES_SELF_EXPANDING_COLLECTION_H_

#include <cassert>

#include "src/data_structures/collection.h"
#include "src/parameters/defined_types.h"

namespace IRL {
/// \brief A modified Collection class
/// that self expands when using the operator[]
/// would force out of memory access.
template <class ObjectType>
class SelfExpandingCollection : public Collection<ObjectType> {
 public:
  /// \brief Default constructor.
  SelfExpandingCollection(void) = default;

  /// \brief This will self-expand to prevent itself from
  /// accessing out of bounds memory.
  ObjectType& operator[](const UnsignedIndex_t a_index);

  /// \brief Const version for access to object in collection.
  const ObjectType& operator[](const UnsignedIndex_t a_index) const;

  /// \brief Default destructor.
  ~SelfExpandingCollection(void) = default;

 private:
};
}  // namespace IRL

#include "src/data_structures/self_expanding_collection.tpp"

#endif  // SRC_DATA_STRUCTURES_SELF_EXPANDING_COLLECTION_H_
