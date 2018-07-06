// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_OPTIMIZATION_BISECTION_TPP_
#define SRC_OPTIMIZATION_BISECTION_TPP_

namespace IRL {
template <class OptimizingClass>
Bisection<OptimizingClass>::Bisection(void) : otype_m(nullptr) {}

template <class OptimizingClass>
void Bisection<OptimizingClass>::solve(OptimizingClass* a_otype,
                                       double a_bracket_0, double a_bracket_1,
                                       const UnsignedIndex_t a_max_iter) {
  otype_m = a_otype;
  this->solve(a_bracket_0, a_bracket_1, a_max_iter);
}

template <class OptimizingClass>
void Bisection<OptimizingClass>::solve(double a_bracket_0, double a_bracket_1,
                                       const UnsignedIndex_t a_max_iter) {
  assert(otype_m != nullptr);

  double error_bracket_0 = this->calculateError(&a_bracket_0);
  if (!otype_m->errorTooHigh(error_bracket_0)) {
    otype_m->setGuess(&a_bracket_0);
    return;
  }
  double error_bracket_1 = this->calculateError(&a_bracket_1);
  if (!otype_m->errorTooHigh(error_bracket_1)) {
    otype_m->setGuess(&a_bracket_1);
    return;
  }
  assert(error_bracket_0 * error_bracket_1 < 0.0);
  double middle = 0.5 * (a_bracket_0 + a_bracket_1);
  double error_middle = this->calculateError(&middle);

  UnsignedIndex_t iteration = {0};
  while (otype_m->errorTooHigh(error_middle)) {
    if (error_bracket_0 * error_middle <= 0.0) {
      a_bracket_1 = middle;
      error_bracket_1 = error_middle;
    } else {
      a_bracket_0 = middle;
      error_bracket_0 = error_middle;
    }
    middle = 0.5 * (a_bracket_0 + a_bracket_1);
    error_middle = this->calculateError(&middle);
    ++iteration;
    if (iteration > a_max_iter) {
      return;
    }
  }
}

template <class OptimizingClass>
double Bisection<OptimizingClass>::calculateError(double* a_guess) {
  otype_m->setGuess(a_guess);
  return otype_m->calculateSignedScalarError();
}

}  // namespace IRL

#endif  // SRC_OPTIMIZATION_BISECTION_TPP_
