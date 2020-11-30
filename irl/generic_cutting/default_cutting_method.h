// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_DEFAULT_CUTTING_METHOD_H_
#define IRL_GENERIC_CUTTING_DEFAULT_CUTTING_METHOD_H_

namespace IRL {

// Cutting methods
struct HalfEdgeCutting {};
struct RecursiveSimplexCutting {};
struct SimplexCutting {};

// Default
using DefaultCuttingMethod = HalfEdgeCutting;
using ReconstructionDefaultCuttingMethod = HalfEdgeCutting;//RecursiveSimplexCutting;
}  // namespace IRL

#endif // IRL_GENERIC_CUTTING_DEFAULT_CUTTING_METHOD_H_
