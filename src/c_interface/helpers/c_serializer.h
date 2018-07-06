// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_C_INTERFACE_HELPERS_C_SERIALIZER_H_
#define SRC_C_INTERFACE_HELPERS_C_SERIALIZER_H_

/// \file c_serializer.h
///
/// These C-style funcions are
/// mapped to functions available in src/serializer.h.
///
/// This file includes functions to handle the serialization and
/// packing of IRL class objects into linear byte-buffers. This also
/// includes the class ByteBuffer, which manages this linear packing
/// and tracks its current buffer location, allowing easy
/// sequential reading that takes place over multiple calls. These functions
/// are mostly planned to be used along with MPI communication routines to
/// send MPI_BYTEs between processors. This means that these functions
/// assume a HOMOGENEOUS ARCHITECTURE, requiring all little-endian or
/// all big-endian representation to be used.

#include "src/c_interface/helpers/c_byte_buffer.h"
#include "src/c_interface/planar_reconstruction/c_separators.h"
#include "src/helpers/serializer.h"
#include "src/parameters/defined_types.h"
#include "src/planar_reconstruction/planar_separator.h"

extern "C" {

void c_serializeAndPack_PlanarSep_ByteBuffer(const c_PlanarSep* a_separator,
                                             c_ByteBuffer* a_container);

void c_unpackAndStore_PlanarSep_ByteBuffer(c_PlanarSep* a_separator,
                                           c_ByteBuffer* a_container);
}

#endif  // SRC_C_INTERFACE_HELPERS_C_SERIALIZER_H_
