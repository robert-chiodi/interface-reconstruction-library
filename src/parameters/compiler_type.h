// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_PARAMETERS_COMPILER_TYPE_H_
#define SRC_PARAMETERS_COMPILER_TYPE_H_

#ifdef __INTEL_COMPILER
#define USING_INTEL_COMPILER
#else
#ifdef __GNUC__
#define USING_GNU_COMPILER
#endif
#endif

#endif  // SRC_PARAMETERS_COMPILER_TYPE_H_
