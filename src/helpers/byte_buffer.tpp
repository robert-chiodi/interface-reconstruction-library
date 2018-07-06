// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_HELPERS_BYTE_BUFFER_TPP_
#define SRC_HELPERS_BYTE_BUFFER_TPP_

namespace IRL {

inline ByteBuffer::ByteBuffer(void) : buffer_m(), buffer_location_m(0) {}

inline LargeOffsetIndex_t ByteBuffer::size(void) const {
  return static_cast<LargeOffsetIndex_t>(buffer_m.size());
}

inline void ByteBuffer::resize(const LargeOffsetIndex_t a_index) {
  buffer_m.resize(a_index);
}

inline Byte_t* ByteBuffer::data(void) { return buffer_m.data(); }

inline const Byte_t* ByteBuffer::data(void) const { return buffer_m.data(); }

inline void ByteBuffer::resetBufferPointer(void) { this->setBufferLocation(0); }

inline void ByteBuffer::reserve(const LargeOffsetIndex_t a_number_to_reserve) {
  assert(a_number_to_reserve < buffer_m.max_size());
  buffer_m.reserve(a_number_to_reserve);
}

template <class ObjectType>
inline void ByteBuffer::pack(const ObjectType* a_list,
                             const LargeOffsetIndex_t a_number_of_elements) {
  const Byte_t* begin = reinterpret_cast<const Byte_t*>(a_list);
  const Byte_t* end = begin + a_number_of_elements * sizeof(ObjectType);
  assert(buffer_location_m <= this->size());
  this->resize(buffer_location_m + a_number_of_elements * sizeof(ObjectType));
  std::copy(begin, end, this->getBufferLocationPointer());
  this->setBufferLocation(static_cast<LargeOffsetIndex_t>(this->size()));
}

template <class ObjectType>
inline void ByteBuffer::unpack(ObjectType* a_list,
                               const LargeOffsetIndex_t a_number_of_elements) {
  Byte_t* object_start = reinterpret_cast<Byte_t*>(a_list);
  const Byte_t* buffer_start = this->getBufferLocationPointer();
  const Byte_t* buffer_end = this->getOffsetBufferLocationPointer(
      a_number_of_elements * sizeof(ObjectType));
  std::copy(buffer_start, buffer_end, object_start);
  this->advanceBufferLocation(static_cast<LargeOffsetIndex_t>(
      a_number_of_elements * sizeof(ObjectType)));
}

inline Byte_t* ByteBuffer::getBufferLocationPointer(void) {
  assert(buffer_location_m < this->size());
  return &buffer_m[buffer_location_m];
}

inline Byte_t* ByteBuffer::getOffsetBufferLocationPointer(
    const LargeOffsetIndex_t a_offset) {
  assert(buffer_location_m + a_offset <= this->size());
  assert(a_offset > 0);
  return &buffer_m[buffer_location_m + a_offset];
}

inline void ByteBuffer::setBufferLocation(const LargeOffsetIndex_t a_index) {
  assert(a_index <= this->size());
  buffer_location_m = a_index;
}

inline void ByteBuffer::advanceBufferLocation(
    const LargeOffsetIndex_t a_increment) {
  buffer_location_m += a_increment;
  assert(buffer_location_m <= this->size());
}

}  // namespace IRL

#endif  // SRC_HELPERS_BYTE_BUFFER_TPP_
