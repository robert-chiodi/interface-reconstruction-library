// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <stdio.h>
#include <sys/stat.h>
#include <string>

#include "examples/paraboloid_advector/diagnostics.h"

void writeOutDiagnostics(const int a_iteration, const int a_revolution,
                         const int a_step, const int a_steps_per_rev) {
  const double percent_complete = 100.0 * static_cast<double>(a_step) /
                                  static_cast<double>(a_steps_per_rev);
  constexpr int bar_length = 20;
  const int number_of_x =
      static_cast<int>(std::ceil(percent_complete)) / (100 / bar_length);

  printf("%14s %4d %10s", "Revolution: ", a_revolution, "|");
  for (int i = 0; i < number_of_x; ++i) {
    printf("x");
  }
  for (int i = number_of_x; i < bar_length; ++i) {
    printf("-");
  }
  printf("|   %6.2f%%", percent_complete);
  printf("\r");
  fflush(stdout);
}

void newlineDiagnostic(void) { printf("\n"); }
