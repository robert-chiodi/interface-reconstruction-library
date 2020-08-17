// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_HELPERS_C_BYTE_BUFFER_H_
#define SRC_C_INTERFACE_HELPERS_C_BYTE_BUFFER_H_

#include "src/helpers/byte_buffer.h"

extern "C" {

struct c_ByteBuffer {
  IRL::ByteBuffer* obj_ptr = nullptr;
};

void c_ByteBuffer_new(c_ByteBuffer* a_self);
void c_ByteBuffer_delete(c_ByteBuffer* a_self);
IRL::LargeOffsetIndex_t c_ByteBuffer_getSize(const c_ByteBuffer* a_self);
void c_ByteBuffer_setSize(c_ByteBuffer* a_self,
                          const IRL::LargeOffsetIndex_t* a_size);
void c_ByteBuffer_resetBufferPointer(c_ByteBuffer* a_self);
IRL::Byte_t* c_ByteBuffer_dataPtr(c_ByteBuffer* a_self);
}

#endif  // SRC_C_INTERFACE_HELPERS_C_BYTE_BUFFER_H_
