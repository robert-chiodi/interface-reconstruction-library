// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_GENERAL_ENCOUNTERED_ID_LIST_H_
#define IRL_GENERIC_CUTTING_GENERAL_ENCOUNTERED_ID_LIST_H_

#include <functional>
#include <utility>

#include "irl/data_structures/unordered_map.h"

#include "irl/data_structures/small_vector.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

class EncounteredIdList {
  using iterator = SmallVector<UnsignedIndex_t, 10>::iterator;
  using const_iterator = SmallVector<UnsignedIndex_t, 10>::const_iterator;

public:
  EncounteredIdList(void) = default;

  UnsignedIndex_t operator[](const UnsignedIndex_t a_index) const {
    assert(a_index < this->size());
    return encountered_list_m[a_index];
  }

  void addId(const UnsignedIndex_t a_id) { encountered_list_m.push_back(a_id); }

  UnsignedIndex_t size(void) const {
    return static_cast<UnsignedIndex_t>(encountered_list_m.size());
  }

  void resize(const UnsignedIndex_t a_size) {
    encountered_list_m.resize(a_size);
  }

  void push_back(const UnsignedIndex_t a_id) {
    encountered_list_m.push_back(a_id);
  }

  void pop_back(void) { encountered_list_m.pop_back(); }

  bool isIdPresent(const UnsignedIndex_t a_id) {
    return std::find(encountered_list_m.begin(), encountered_list_m.end(),
                     a_id) != encountered_list_m.end();
  }

  iterator begin(void) { return encountered_list_m.begin(); }
  const_iterator begin(void) const { return this->cbegin(); }
  const_iterator cbegin(void) const { return encountered_list_m.cbegin(); }
  iterator end(void) { return encountered_list_m.end(); }
  const_iterator end(void) const { return this->cend(); }
  const_iterator cend(void) const { return encountered_list_m.cend(); }

  ~EncounteredIdList(void) = default;

private:
  SmallVector<UnsignedIndex_t, 10> encountered_list_m;
};

// class EncounteredIdList {
//   static constexpr std::size_t initial_capacity = 10;

//  public:
//   EncounteredIdList(void) : valid_size_m(0), encountered_list_m() {
//     // encountered_list_m.reserve(initial_capacity);
//   }

//   void addId(const UnsignedIndex_t a_id) {
//     encountered_list_m[a_id] = this->size();
//     ++valid_size_m;
//   }

//   UnsignedIndex_t size(void) const { return valid_size_m; }

//   void resize(const UnsignedIndex_t a_size) { valid_size_m = a_size; }

//   bool isIdPresent(const UnsignedIndex_t a_id) {
//     return encountered_list_m.find(a_id) != encountered_list_m.end()
//                ? encountered_list_m.at(a_id) < valid_size_m
//                : false;
//   }

//   ~EncounteredIdList(void) = default;

//  private:
//   UnsignedIndex_t valid_size_m;
//   IRL::unordered_map<UnsignedIndex_t, UnsignedIndex_t> encountered_list_m;
// };

} // namespace IRL

#endif // SRC_GENERIC_CUTTING_GENERAL_ENCOUNTERED_ID_LIST_H_
