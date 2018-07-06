// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_LISTER_H_
#define SRC_DATA_STRUCTURES_LISTER_H_

#include "src/data_structures/collection.h"
#include "src/parameters/defined_types.h"

namespace IRL {

/// \brief This class is a Collection that
/// has operator += add the other Lister or
/// object to the back of its own collection.
template <class ObjectType>
class Lister : public Collection<ObjectType> {
 public:
  /// \brief Default constructor.
  Lister(void) = default;

  /// \brief The operator += will be used to append
  /// `a_rhs` to the current collection in Lister.
  Lister& operator+=(const Lister& a_rhs);

  /// \brief The operator += will be used to push_back
  /// the object `a_rhs` in the current collection in Lister.
  Lister& operator+=(const ObjectType& a_rhs);

  /// \brief Default destructor.
  ~Lister(void) = default;

 private:
};

}  // namespace IRL

#include "src/data_structures/lister.tpp"

#endif  // SRC_DATA_STRUCTURES_LISTER_H_
