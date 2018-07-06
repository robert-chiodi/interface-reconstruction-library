// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_OBJECT_ALLOCATION_SERVER_H_
#define SRC_DATA_STRUCTURES_OBJECT_ALLOCATION_SERVER_H_

#include "src/data_structures/block_object_allocation.h"
#include "src/parameters/defined_types.h"

namespace IRL {

template <class ObjectType>
class ObjectAllocationServer {
 public:
  ObjectAllocationServer(void) = delete;

  explicit ObjectAllocationServer(
      const LargeOffsetIndex_t a_number_to_allocate);

  ObjectType* getNewObject(void);

  explicit ObjectAllocationServer(const ObjectAllocationServer& other) = delete;

  explicit ObjectAllocationServer(ObjectAllocationServer&& other) noexcept =
      delete;

  ObjectAllocationServer& operator=(const ObjectAllocationServer& other) =
      delete;

  ObjectAllocationServer& operator=(ObjectAllocationServer&& other) noexcept =
      delete;

  ~ObjectAllocationServer(void) = default;

 private:
  ObjectType* serveNextAvailableObject(void);

  LargeOffsetIndex_t getIndexForNextObjectToServe(void);

  void incrementNumberOfServedObjects(void);

  LargeOffsetIndex_t indexOfNextObjectToServe(void) const;

  BlockObjectAllocation<ObjectType> object_allocation_m;
  LargeOffsetIndex_t number_of_objects_served_m;
};
}  // namespace IRL

#include "src/data_structures/object_allocation_server.tpp"

#endif  // SRC_DATA_STRUCTURES_OBJECT_ALLOCATION_SERVER_H_
