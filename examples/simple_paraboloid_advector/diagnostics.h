// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_SIMPLE_PARABOLOID_ADVECTOR_DIAGNOSTICS_H_
#define EXAMPLES_SIMPLE_PARABOLOID_ADVECTOR_DIAGNOSTICS_H_

#include "examples/simple_paraboloid_advector/basic_mesh.h"
#include "examples/simple_paraboloid_advector/data.h"

void writeOutDiagnostics(const int a_iteration, const int a_revolutions,
                         const int a_step, const int a_steps_per_rev);

void newlineDiagnostic(void);

#endif  // EXAMPLES_SIMPLE_PARABOLOID_ADVECTOR_DIAGNOSTICS_H_
