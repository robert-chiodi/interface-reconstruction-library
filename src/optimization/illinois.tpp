// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_OPTIMIZATION_ILLINOIS_TPP_
#define SRC_OPTIMIZATION_ILLINOIS_TPP_

namespace IRL {

template <class OptimizingClass>
Illinois<OptimizingClass>::Illinois(void) : otype_m(nullptr) {}

template <class OptimizingClass>
void Illinois<OptimizingClass>::solve(OptimizingClass* a_otype,
                                      double a_bracket_0, double a_bracket_1) {
  otype_m = a_otype;
  this->solve(a_bracket_0, a_bracket_1);
}

template <class OptimizingClass>
void Illinois<OptimizingClass>::solve(double a_bracket_0, double a_bracket_1) {
  assert(otype_m != nullptr);

  double error_bracket_0 = this->calculateError(&a_bracket_0);
  double error_bracket_1 = this->calculateError(&a_bracket_1);
  assert(error_bracket_0 * error_bracket_1 < 0.0);

  int side = 0;
  double middle, error_middle = DBL_MAX;
  while (otype_m->errorTooHigh(error_middle)) {
    middle = (error_bracket_0 * a_bracket_1 - error_bracket_1 * a_bracket_0) /
             (error_bracket_0 - error_bracket_1);
    error_middle = this->calculateError(&middle);

    if (error_bracket_0 * error_middle <= 0.0) {
      a_bracket_1 = middle;
      error_bracket_1 = error_middle;
      if (side == -1) {
        error_bracket_0 *= 0.5;
      }
      side = -1;
    } else {
      a_bracket_0 = middle;
      error_bracket_0 = error_middle;
      if (side == 1) {
        error_bracket_1 *= 0.5;
      }
      side = 1;
    }
  }
}

template <class OptimizingClass>
double Illinois<OptimizingClass>::calculateError(double* a_guess) {
  otype_m->setGuess(a_guess);
  return otype_m->calculateSignedScalarError();
}

}  // namespace IRL

#endif  // SRC_OPTIMIZATION_ILLINOIS_TPP_
