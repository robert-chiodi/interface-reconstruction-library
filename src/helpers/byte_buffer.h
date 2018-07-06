// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_HELPERS_BYTE_BUFFER_H_
#define SRC_HELPERS_BYTE_BUFFER_H_

#include <algorithm>
#include <cassert>
#include <vector>

#include "src/parameters/defined_types.h"

namespace IRL {

class ByteBuffer {
 public:
  /// \brief Default constructor.
  ByteBuffer(void);

  /// \brief Return const size of the buffer.
  LargeOffsetIndex_t size(void) const;

  /// \brief Resize the underlying buffer array.
  void resize(const LargeOffsetIndex_t a_index);

  /// \brief Return object pointer to at start of buffer_m.
  Byte_t* data(void);

  /// \brief Return const object pointer to at start of buffer_m.
  const Byte_t* data(void) const;

  /// \brief Resets the buffer location index to 0 (the start).
  void resetBufferPointer(void);

  /// \brief Allow reservation of space to gain better performance if
  /// size is known a priori.
  void reserve(const LargeOffsetIndex_t a_number_to_reserve);

  /// \brief Pack object (or list of objects) onto the
  /// end of the buffer.
  template <class ObjectType>
  void pack(const ObjectType* a_list,
            const LargeOffsetIndex_t a_number_of_elements);

  /// \brief Unpack bytes in list, starting at current buffer location
  /// into objects pointed to by `a_list`.
  template <class ObjectType>
  void unpack(ObjectType* a_list,
              const LargeOffsetIndex_t a_number_of_elements);

  /// \brief Default destructor.
  ~ByteBuffer(void) = default;

 private:
  /// \brief Get pointer to underlying data at current buffer location.
  Byte_t* getBufferLocationPointer(void);

  /// \brief Get pointer to underlying data at offset from current buffer
  /// location.
  Byte_t* getOffsetBufferLocationPointer(const LargeOffsetIndex_t a_offset);

  /// \brief Set buffer location index to `a_index`.
  void setBufferLocation(const LargeOffsetIndex_t a_index);

  /// \brief Update buffer location.
  void advanceBufferLocation(const LargeOffsetIndex_t a_increment);

  /// \brief Vector of bytes.
  std::vector<Byte_t> buffer_m;
  /// \brief Current buffer location.
  // Done as index since expanding vector might invalidate a pointer.
  LargeOffsetIndex_t buffer_location_m;
};

}  // namespace IRL

#include "src/helpers/byte_buffer.tpp"

#endif  // SRC_HELPERS_BYTE_BUFFER_H_
