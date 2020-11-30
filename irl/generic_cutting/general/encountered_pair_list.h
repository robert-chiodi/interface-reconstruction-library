// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_GENERAL_ENCOUNTERED_PAIR_LIST_H_
#define IRL_GENERIC_CUTTING_GENERAL_ENCOUNTERED_PAIR_LIST_H_

#include <functional>
#include <iostream>
#include <utility>

#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/unordered_map.h"

#include "irl/parameters/defined_types.h"

namespace IRL {

template <class ContainedType>
class TmpBadName {
  static constexpr UnsignedIndex_t max_val = static_cast<UnsignedIndex_t>(-1);
  struct Starts {
    Starts(void) noexcept
        : key_m(max_val),
          size_m(max_val),
          start_index_m(max_val),
          end_index_m(max_val) {}

    Starts(const UnsignedIndex_t a_key, const UnsignedIndex_t a_size,
           const UnsignedIndex_t a_start_index,
           const UnsignedIndex_t a_end_index) noexcept
        : key_m(a_key),
          size_m(a_size),
          start_index_m(a_start_index),
          end_index_m(a_end_index) {}
    UnsignedIndex_t key_m, size_m, start_index_m, end_index_m;
  };
  struct LinkedListMember {
    LinkedListMember(void) noexcept : value_m(), next_m(max_val) {}

    LinkedListMember(const ContainedType &a_value,
                     const UnsignedIndex_t a_next) noexcept
        : value_m(a_value), next_m(a_next) {}
    ContainedType value_m;
    UnsignedIndex_t next_m;
  };

 public:
  struct State {
    State(const std::vector<Starts> &a_start_vector,
          const UnsignedIndex_t a_linked_list_size) noexcept
        : starts_vector_m(a_start_vector),
          linked_list_storage_size_m(a_linked_list_size) {}

    std::vector<Starts> starts_vector_m;
    UnsignedIndex_t linked_list_storage_size_m;
  };

 public:
  void addMember(const UnsignedIndex_t a_key, const ContainedType &a_value) {
    auto location = this->findStartingLocation(a_key);
    UnsignedIndex_t new_value_location =
        static_cast<UnsignedIndex_t>(linked_list_storage_m.size());
    linked_list_storage_m.emplace_back(a_value,
                                       static_cast<UnsignedIndex_t>(-1));
    // Location not present and is largest key value
    if (location == starting_locations_m.end()) {
      starting_locations_m.emplace_back(a_key, 1, new_value_location,
                                        new_value_location);

      // Key not present but needs to be inserted to maintain ascending
      // order
    } else if ((*location).key_m != a_key) {
      starting_locations_m.emplace(location, a_key, 1, new_value_location,
                                   new_value_location);

      // Key present, need to add item to linked list
    } else {
      linked_list_storage_m[(*location).end_index_m].next_m =
          new_value_location;
      (*location).end_index_m = new_value_location;
      ++(*location).size_m;
    }
  }

  bool isPairPresent(const UnsignedIndex_t a_key,
                     const ContainedType &a_value) const {
    auto location = this->findStartingLocation(a_key);
    if (location == starting_locations_m.end() || (*location).key_m != a_key) {
      return false;
    }

    return this->isPresentInLinkedList((*location).start_index_m,
                                       (*location).size_m, a_value);
  }

  bool addIfNotPresent(const UnsignedIndex_t a_key,
                       const ContainedType &a_value) {
    auto location = this->findStartingLocation(a_key);
    if (location != starting_locations_m.end()) {
      if ((*location).key_m == a_key) {
        if (this->isPresentInLinkedList((*location).start_index_m,
                                        (*location).size_m, a_value)) {
          return true;
        }
      }
    }

    // Otherwise this is a new pair, combine
    UnsignedIndex_t new_value_location =
        static_cast<UnsignedIndex_t>(linked_list_storage_m.size());
    linked_list_storage_m.emplace_back(a_value,
                                       static_cast<UnsignedIndex_t>(-1));
    // Location not present and is largest key value
    if (location == starting_locations_m.end()) {
      starting_locations_m.emplace_back(a_key, 1, new_value_location,
                                        new_value_location);

      // Key not present but needs to be inserted to maintain ascending
      // order
    } else if ((*location).key_m != a_key) {
      starting_locations_m.emplace(location, a_key, 1, new_value_location,
                                   new_value_location);

      // Key present, need to add item to linked list
    } else {
      linked_list_storage_m[(*location).end_index_m].next_m =
          new_value_location;
      (*location).end_index_m = new_value_location;
      ++(*location).size_m;
    }
    return false;
  }

  // Export a vector that can be used to efficiently reset this to a
  // previous state.
  State saveState(void) const {
    return {starting_locations_m,
            static_cast<UnsignedIndex_t>(linked_list_storage_m.size())};
  }

  void resetToState(const State &a_previous_state) {
    linked_list_storage_m.resize(a_previous_state.linked_list_storage_size_m);
    starting_locations_m.resize(a_previous_state.starts_vector_m.size());
    std::copy(a_previous_state.starts_vector_m.begin(),
              a_previous_state.starts_vector_m.end(),
              starting_locations_m.begin());
  }

  void printOut(void) const {
    for (const auto &member : starting_locations_m) {
      std::cout << "ID of " << member.key_m << " connected to ";
      UnsignedIndex_t index = member.start_index_m;
      for (UnsignedIndex_t n = 0; n < member.size_m; ++n) {
        std::cout << linked_list_storage_m[index].value_m << " ";
        index = linked_list_storage_m[index].next_m;
      }
      std::cout << '\n' << std::endl;
    }
  }

 private:
  // starting_locations_m is assumed to be sorted in ascending order.
  typename std::vector<Starts>::iterator findStartingLocation(
      const UnsignedIndex_t a_key) {
    if (!starting_locations_m.empty() &&
        a_key > starting_locations_m.back().key_m) {
      return starting_locations_m.end();
    } else if (!starting_locations_m.empty() &&
               a_key <= starting_locations_m.front().key_m) {
      return starting_locations_m.begin();
    } else {
      UnsignedIndex_t lower = {0};
      UnsignedIndex_t upper =
          static_cast<UnsignedIndex_t>(starting_locations_m.size());
      while (upper - lower > 1) {
        UnsignedIndex_t middle = (upper + lower) / 2;
        if (starting_locations_m[middle].key_m == a_key) {
          return starting_locations_m.begin() + middle;
        } else if (starting_locations_m[middle].key_m < a_key) {
          lower = middle;
        } else {
          upper = middle;
        }
      }
      return starting_locations_m.begin() + upper;
    }
  }

  typename std::vector<Starts>::const_iterator findStartingLocation(
      const UnsignedIndex_t a_key) const {
    if (!starting_locations_m.empty() &&
        a_key > starting_locations_m.back().key_m) {
      return starting_locations_m.end();
    } else if (!starting_locations_m.empty() &&
               a_key <= starting_locations_m.front().key_m) {
      return starting_locations_m.begin();
    } else {
      UnsignedIndex_t lower = {0};
      UnsignedIndex_t upper =
          static_cast<UnsignedIndex_t>(starting_locations_m.size());
      while (upper - lower > 1) {
        UnsignedIndex_t middle = (upper + lower) / 2;
        if (starting_locations_m[middle].key_m == a_key) {
          return starting_locations_m.begin() + middle;
        } else if (starting_locations_m[middle].key_m < a_key) {
          lower = middle;
        } else {
          upper = middle;
        }
      }
      return starting_locations_m.begin() + upper;
    }
  }

  bool isPresentInLinkedList(const UnsignedIndex_t a_start,
                             const UnsignedIndex_t a_size,
                             const ContainedType &a_value) const {
    UnsignedIndex_t index = a_start;
    for (UnsignedIndex_t n = 0; n < a_size; ++n) {
      assert(index < linked_list_storage_m.size());
      if (linked_list_storage_m[index].value_m == a_value) {
        return true;
      }
      index = linked_list_storage_m[index].next_m;
    }
    return false;
  }

  std::vector<Starts> starting_locations_m;
  std::vector<LinkedListMember> linked_list_storage_m;
};

template <class ContainedType>
class TmpBadNameNotSorted {
  static constexpr UnsignedIndex_t max_val = static_cast<UnsignedIndex_t>(-1);
  struct Starts {
    Starts(void) noexcept
        : key_m(max_val),
          size_m(max_val),
          start_index_m(max_val),
          end_index_m(max_val) {}

    Starts(const UnsignedIndex_t a_key, const UnsignedIndex_t a_size,
           const UnsignedIndex_t a_start_index,
           const UnsignedIndex_t a_end_index) noexcept
        : key_m(a_key),
          size_m(a_size),
          start_index_m(a_start_index),
          end_index_m(a_end_index) {}
    UnsignedIndex_t key_m, size_m, start_index_m, end_index_m;
  };
  struct LinkedListMember {
    LinkedListMember(void) noexcept : value_m(), next_m(max_val) {}

    LinkedListMember(const ContainedType &a_value,
                     const UnsignedIndex_t a_next) noexcept
        : value_m(a_value), next_m(a_next) {}
    ContainedType value_m;
    UnsignedIndex_t next_m;
  };

 public:
  using State = std::vector<UnsignedIndex_t>;

 public:
  void addMember(const UnsignedIndex_t a_key, const ContainedType &a_value) {
    auto location = starting_locations_m.begin();
    while (location != starting_locations_m.end()) {
      if ((*location).key_m == a_key) {
        break;
      }
      ++location;
    }
    UnsignedIndex_t new_value_location =
        static_cast<UnsignedIndex_t>(linked_list_storage_m.size());
    linked_list_storage_m.emplace_back(a_value,
                                       static_cast<UnsignedIndex_t>(-1));
    // Location not present and is largest key value
    if (location == starting_locations_m.end()) {
      starting_locations_m.emplace_back(a_key, 1, new_value_location,
                                        new_value_location);

      // Key present, need to add item to linked list
    } else {
      linked_list_storage_m[(*location).end_index_m].next_m =
          new_value_location;
      (*location).end_index_m = new_value_location;
      ++(*location).size_m;
    }
  }

  bool isPairPresent(const UnsignedIndex_t a_key,
                     const ContainedType &a_value) const {
    auto location = starting_locations_m.begin();
    while (location != starting_locations_m.end()) {
      if ((*location).key_m == a_key) {
        break;
      }
      ++location;
    }
    if (location != starting_locations_m.end()) {
      if (this->isPresentInLinkedList((*location).start_index_m,
                                      (*location).size_m, a_value)) {
        return true;
      }
    }
    return false;
  }

  bool addIfNotPresent(const UnsignedIndex_t a_key,
                       const ContainedType &a_value) {
    auto location = starting_locations_m.begin();
    while (location != starting_locations_m.end()) {
      if ((*location).key_m == a_key) {
        break;
      }
      ++location;
    }
    if (location != starting_locations_m.end()) {
      if (this->isPresentInLinkedList((*location).start_index_m,
                                      (*location).size_m, a_value)) {
        return true;
      }
    }

    // Otherwise this is a new pair, combine
    UnsignedIndex_t new_value_location =
        static_cast<UnsignedIndex_t>(linked_list_storage_m.size());
    linked_list_storage_m.emplace_back(a_value,
                                       static_cast<UnsignedIndex_t>(-1));
    // Location not present
    if (location == starting_locations_m.end()) {
      starting_locations_m.emplace_back(a_key, 1, new_value_location,
                                        new_value_location);
      // Key present, need to add item to linked list
    } else {
      linked_list_storage_m[(*location).end_index_m].next_m =
          new_value_location;
      (*location).end_index_m = new_value_location;
      ++(*location).size_m;
    }
    return false;
  }

  // Export a vector that can be used to efficiently reset this to a
  // previous state.
  State saveState(void) const {
    std::vector<UnsignedIndex_t> state(2 + starting_locations_m.size());
    state[0] = static_cast<UnsignedIndex_t>(starting_locations_m.size());
    state[1] = static_cast<UnsignedIndex_t>(linked_list_storage_m.size());
    auto state_iterator = state.begin() + 2;
    for (const auto &member : starting_locations_m) {
      (*state_iterator) = member.size_m;
      ++state_iterator;
    }
    return state;
  }

  void resetToState(const State &a_previous_state) {
    starting_locations_m.resize(a_previous_state[0]);
    linked_list_storage_m.resize(a_previous_state[1]);
    auto state_iterator = a_previous_state.begin() + 2;
    for (auto &member : starting_locations_m) {
      member.size_m = (*state_iterator);
      ++state_iterator;
    }

    for (auto &member : starting_locations_m) {
      auto index = member.start_index_m;
      for (UnsignedIndex_t n = 1; n < member.size_m; ++n) {
        index = linked_list_storage_m[index].next_m;
      }
      assert(index < linked_list_storage_m.size());
      member.end_index_m = index;
    }
  }

  void printOut(void) const {
    for (const auto &member : starting_locations_m) {
      std::cout << "ID of " << member.key_m << " connected to ";
      UnsignedIndex_t index = member.start_index_m;
      for (UnsignedIndex_t n = 0; n < member.size_m; ++n) {
        std::cout << linked_list_storage_m[index].value_m << " ";
        index = linked_list_storage_m[index].next_m;
      }
      std::cout << '\n' << std::endl;
    }
  }

 private:
  bool isPresentInLinkedList(const UnsignedIndex_t a_start,
                             const UnsignedIndex_t a_size,
                             const ContainedType &a_value) const {
    UnsignedIndex_t index = a_start;
    for (UnsignedIndex_t n = 0; n < a_size; ++n) {
      assert(index < linked_list_storage_m.size());
      if (linked_list_storage_m[index].value_m == a_value) {
        return true;
      }
      index = linked_list_storage_m[index].next_m;
    }
    return false;
  }

  std::vector<Starts> starting_locations_m;
  std::vector<LinkedListMember> linked_list_storage_m;
};

class EncounteredPairList {
  using PairStorage = TmpBadNameNotSorted<UnsignedIndex_t>;

 public:
  EncounteredPairList(void) : encountered_list_m() {
    // encountered_list_m.reserve(initial_capacity);
  }

  void addPair(const UnsignedIndex_t a_id_0, const UnsignedIndex_t a_id_1) {
    assert(a_id_0 != a_id_1);
    const auto lower_id = a_id_0 < a_id_1 ? a_id_0 : a_id_1;
    const auto upper_id = a_id_0 > a_id_1 ? a_id_0 : a_id_1;
    encountered_list_m.addMember(lower_id, upper_id);
  }

  bool isPairPresent(const UnsignedIndex_t a_id_0,
                     const UnsignedIndex_t a_id_1) const {
    const auto lower_id = a_id_0 < a_id_1 ? a_id_0 : a_id_1;
    const auto upper_id = a_id_0 > a_id_1 ? a_id_0 : a_id_1;
    return encountered_list_m.isPairPresent(lower_id, upper_id);
  }

  bool addIfNotPresent(const UnsignedIndex_t a_id_0,
                       const UnsignedIndex_t a_id_1) {
    const auto lower_id = a_id_0 < a_id_1 ? a_id_0 : a_id_1;
    const auto upper_id = a_id_0 > a_id_1 ? a_id_0 : a_id_1;
    return encountered_list_m.addIfNotPresent(lower_id, upper_id);
  }

  auto saveState(void) const { return encountered_list_m.saveState(); }

  void resetToState(const typename PairStorage::State &a_previous_state) {
    encountered_list_m.resetToState(a_previous_state);
  }

  void printOut(void) const { encountered_list_m.printOut(); }

  ~EncounteredPairList(void) = default;

 private:
  PairStorage encountered_list_m;
};

// class EncounteredPairList {
//   static constexpr std::size_t initial_capacity = 10;

// public:
//   EncounteredPairList(void) : encountered_list_m() {
//     // encountered_list_m.reserve(initial_capacity);
//   }

//   void addPair(const UnsignedIndex_t a_id_0, const UnsignedIndex_t a_id_1)
//   {
//     assert(a_id_0 != a_id_1);
//     const auto lower_id = a_id_0 < a_id_1 ? a_id_0 : a_id_1;
//     const auto upper_id = a_id_0 > a_id_1 ? a_id_0 : a_id_1;
//     auto loc = encountered_list_m.find(lower_id);
//     if (loc == encountered_list_m.end()) {
//       auto insert_loc = encountered_list_m.insert(
//           {lower_id, SmallVector<UnsignedIndex_t, 6>()});
//       assert(insert_loc.second);
//       ((*(insert_loc.first)).second).push_back(upper_id);
//     } else {
//       (*loc).second.push_back(upper_id);
//     }
//   }

//   bool isPairPresent(const UnsignedIndex_t a_id_0,
//                      const UnsignedIndex_t a_id_1) const {
//     const auto lower_id = a_id_0 < a_id_1 ? a_id_0 : a_id_1;
//     const auto upper_id = a_id_0 > a_id_1 ? a_id_0 : a_id_1;
//     auto loc = encountered_list_m.find(lower_id);
//     if (loc == encountered_list_m.end()) {
//       return false;
//     }
//     const auto &vec = (*loc).second;
//     auto vec_loc = std::find(vec.begin(), vec.end(), upper_id);
//     return vec_loc != vec.end();
//   }

//   void printOut(void) const {
//     for (const auto &elem : encountered_list_m) {
//       std::cout << elem.first << ":";
//       for (const auto &vec_elem : elem.second) {
//         std::cout << " " << vec_elem;
//       }
//       std::cout << std::endl;
//     }
//   }

//   ~EncounteredPairList(void) = default;

// private:
//   IRL::unordered_map<UnsignedIndex_t, SmallVector<UnsignedIndex_t, 6>>
//       encountered_list_m;
// };

}  // namespace IRL

#endif // IRL_GENERIC_CUTTING_GENERAL_ENCOUNTERED_PAIR_LIST_H_
