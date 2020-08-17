// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/c_interface/helpers/c_serializer.h"

#include <cassert>

extern "C" {

void c_serializeAndPack_PlanarSep_ByteBuffer(const c_PlanarSep* a_separator,
                                             c_ByteBuffer* a_container) {
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_container != nullptr);
  assert(a_container->obj_ptr != nullptr);
  IRL::serializeAndPack(*a_separator->obj_ptr, a_container->obj_ptr);
}

void c_unpackAndStore_PlanarSep_ByteBuffer(c_PlanarSep* a_separator,
                                           c_ByteBuffer* a_container) {
  assert(a_separator != nullptr);
  assert(a_separator->obj_ptr != nullptr);
  assert(a_container != nullptr);
  assert(a_container->obj_ptr != nullptr);
  IRL::unpackAndStore(a_separator->obj_ptr, a_container->obj_ptr);
}

}  // end extern C
