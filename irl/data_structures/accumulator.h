// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_DATA_STRUCTURES_ACCUMULATOR_H_
#define IRL_DATA_STRUCTURES_ACCUMULATOR_H_

#include "irl/data_structures/self_expanding_collection.h"

namespace IRL {
/// \brief This class is a SelfExpandingCollection that
/// has operator += call the += operator for each
/// member in its collection.
template <class ObjectType>
class Accumulator : public SelfExpandingCollection<ObjectType> {
 public:
  /// \brief Default constructor.
  Accumulator(void) = default;

  /// \brief The operator += will be used to add another
  /// Accumulator, where the entry in each index is summed
  /// between the two.
  Accumulator& operator+=(const Accumulator& a_rhs);

  /// \brief Default destructor.
  ~Accumulator(void) = default;

 private:
};
}  // namespace IRL

#include "irl/data_structures/accumulator.tpp"

#endif // IRL_DATA_STRUCTURES_ACCUMULATOR_H_
