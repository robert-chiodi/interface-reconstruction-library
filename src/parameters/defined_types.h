// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PARAMETERS_DEFINED_TYPES_H_
#define SRC_PARAMETERS_DEFINED_TYPES_H_

#include <cstddef>
#include <cstdint>

namespace IRL {
using UnsignedIndex_t = uint32_t;
using LargeOffsetIndex_t = std::size_t;
using LookupIndex_t = uint8_t;
using Byte_t = unsigned char;
}  // namespace IRL

#endif  // SRC_DEFINED_TYPES_H_
