// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/helpers/c_byte_buffer.h"

#include <iostream>

extern "C" {

void c_ByteBuffer_new(c_ByteBuffer* a_self) {
  assert(a_self->obj_ptr == nullptr);
  a_self->obj_ptr = new IRL::ByteBuffer;
}

void c_ByteBuffer_delete(c_ByteBuffer* a_self) {
  delete a_self->obj_ptr;
  a_self->obj_ptr = nullptr;
}

IRL::LargeOffsetIndex_t c_ByteBuffer_getSize(const c_ByteBuffer* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return a_self->obj_ptr->size();
}

void c_ByteBuffer_setSize(c_ByteBuffer* a_self,
                          const IRL::LargeOffsetIndex_t* a_size) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->resize(*a_size);
}

void c_ByteBuffer_resetBufferPointer(c_ByteBuffer* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  a_self->obj_ptr->resetBufferPointer();
}

IRL::Byte_t* c_ByteBuffer_dataPtr(c_ByteBuffer* a_self) {
  assert(a_self != nullptr);
  assert(a_self->obj_ptr != nullptr);
  return a_self->obj_ptr->data();
}
}
