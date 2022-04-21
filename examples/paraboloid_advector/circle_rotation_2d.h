// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_PARABOLOID_ADVECTOR_CIRCLE_ROTATION_2D_H_
#define EXAMPLES_PARABOLOID_ADVECTOR_CIRCLE_ROTATION_2D_H_

#include "irl/paraboloid_reconstruction/paraboloid.h"

#include "examples/paraboloid_advector/basic_mesh.h"
#include "examples/paraboloid_advector/data.h"

struct CircleRotation2D {
  static BasicMesh setMesh(void);

  static void initialize(Data<double>* a_U, Data<double>* a_V,
                         Data<double>* a_W, Data<IRL::Paraboloid>* a_interface);

  static void setVelocity(const double a_time, Data<double>* a_U,
                          Data<double>* a_V, Data<double>* a_W);
};

#endif  // EXAMPLES_PARABOLOID_ADVECTOR_CIRCLE_ROTATION_2D_H_
