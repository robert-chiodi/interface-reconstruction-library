// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2023 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef TOOLS_SRC_SURFACE_IO_SURFACE_IO_H_
#define TOOLS_SRC_SURFACE_IO_SURFACE_IO_H_

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "irl/paraboloid_reconstruction/parametrized_surface.h"
#include "irl/parameters/defined_types.h"
#include "irl/surface_mesher/triangulated_surface.h"

int main(int argc, char* argv[]);

void simpleErrorHandler(const std::string& error_message);

void simpleErrorHandler(const std::string& error_message) {
  std::cout << error_message << std::endl;
  exit(1);
}

#endif  // TOOLS_SRC_SURFACE_IO_SURFACE_IO_H_
