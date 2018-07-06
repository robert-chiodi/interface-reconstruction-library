// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_DATA_STRUCTURES_OBJECT_ALLOCATION_SERVER_TPP_
#define SRC_DATA_STRUCTURES_OBJECT_ALLOCATION_SERVER_TPP_

namespace IRL {

template <class ObjectType>
ObjectAllocationServer<ObjectType>::ObjectAllocationServer(
    const LargeOffsetIndex_t a_number_to_allocate)
    : object_allocation_m{a_number_to_allocate},
      number_of_objects_served_m{0} {}

template <class ObjectType>
ObjectType* ObjectAllocationServer<ObjectType>::getNewObject(void) {
  return this->serveNextAvailableObject();
}

template <class ObjectType>
ObjectType* ObjectAllocationServer<ObjectType>::serveNextAvailableObject(void) {
  const LargeOffsetIndex_t index_of_next_available_object =
      this->getIndexForNextObjectToServe();
  return object_allocation_m.getObjectAtIndex(index_of_next_available_object);
}

template <class ObjectType>
LargeOffsetIndex_t
ObjectAllocationServer<ObjectType>::getIndexForNextObjectToServe(void) {
  this->incrementNumberOfServedObjects();
  return indexOfNextObjectToServe();
}

template <class ObjectType>
void ObjectAllocationServer<ObjectType>::incrementNumberOfServedObjects(void) {
  ++number_of_objects_served_m;
}

template <class ObjectType>
LargeOffsetIndex_t ObjectAllocationServer<ObjectType>::indexOfNextObjectToServe(
    void) const {
  return number_of_objects_served_m - 1;
}

}  // namespace IRL

#endif  // SRC_DATA_STRUCTURES_OBJECT_ALLOCATION_SERVER_TPP_
