// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_HELPERS_SFINAE_BOILER_PLATE_H_
#define SRC_HELPERS_SFINAE_BOILER_PLATE_H_

#include <type_traits>

namespace IRL {

template <bool Cond, typename T = void>
using enable_if_t = typename std::enable_if<Cond, T>::type;

}  // namespace IRL

#endif  // SRC_HELPERS_SFINAE_BOILER_PLATE_H_
