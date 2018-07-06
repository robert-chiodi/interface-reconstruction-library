// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_OPTIMIZATION_BRENTS_METHOD_TPP_
#define SRC_OPTIMIZATION_BRENTS_METHOD_TPP_

namespace IRL {

template <class OptimizingClass>
BrentsMethod<OptimizingClass>::BrentsMethod(void) : otype_m(nullptr) {}

template <class OptimizingClass>
void BrentsMethod<OptimizingClass>::solve(OptimizingClass* a_otype,
                                          double a_bracket_0,
                                          double a_bracket_1,
                                          const double a_delta) {
  otype_m = a_otype;
  this->solve(a_bracket_0, a_bracket_1, a_delta);
}

template <class OptimizingClass>
void BrentsMethod<OptimizingClass>::solve(double a_bracket_0,
                                          double a_bracket_1,
                                          const double a_delta) {
  assert(otype_m != nullptr);

  assert(a_delta > 0.0);
  double error_bracket_0 = {calculateError(&a_bracket_0)};
  double error_bracket_1 = {calculateError(&a_bracket_1)};
  assert(error_bracket_0 * error_bracket_1 < 0.0);
  if (std::fabs(error_bracket_0) < std::fabs(error_bracket_1)) {
    std::swap(a_bracket_0, a_bracket_1);
    std::swap(error_bracket_0, error_bracket_1);
  }
  double c = {a_bracket_0};
  double error_c = {error_bracket_0};
  double d = {0.0};
  double s = {0.0};
  double error_s = {error_bracket_1};
  bool mflag = {true};
  double error = std::fabs(error_s);
  while (otype_m->errorTooHigh(error)) {
    if (error_bracket_0 != error_c && error_bracket_1 != error_c) {
      s = a_bracket_0 * error_bracket_1 * error_c /
              ((error_bracket_0 - error_bracket_1) *
               (error_bracket_0 - error_c)) +
          a_bracket_1 * error_bracket_0 * error_c /
              ((error_bracket_1 - error_bracket_0) *
               (error_bracket_1 - error_c)) +
          c * error_bracket_0 * error_bracket_1 /
              ((error_c - error_bracket_0) * (error_c - error_bracket_1));
    } else {
      s = a_bracket_1 - error_bracket_1 * (a_bracket_1 - a_bracket_0) /
                            (error_bracket_1 - error_bracket_0);
    }
    if (s < 0.25 * (3.0 * a_bracket_0 + a_bracket_1) || s > a_bracket_1 ||
        (mflag &&
         std::fabs(s - a_bracket_1) >= 0.5 * std::fabs(a_bracket_1 - c)) ||
        (!mflag && std::fabs(s - a_bracket_1) >= 0.5 * std::fabs(c - d)) ||
        (mflag && std::fabs(a_bracket_1 - c) < a_delta) ||
        (!mflag && std::fabs(c - d) < a_delta)) {
      s = 0.5 * (a_bracket_0 + a_bracket_1);
      mflag = true;
    } else {
      mflag = false;
    }
    error_s = this->calculateError(&s);
    d = c;
    c = a_bracket_1;
    error_c = error_bracket_1;
    if (error_bracket_0 * error_s < 0.0) {
      a_bracket_1 = s;
      error_bracket_1 = error_s;
    } else {
      a_bracket_0 = s;
      error_bracket_0 = error_s;
    }
    if (std::fabs(error_bracket_0) < std::fabs(error_bracket_1)) {
      std::swap(a_bracket_0, a_bracket_1);
      std::swap(error_bracket_0, error_bracket_1);
    }
    error = {std::min(std::fabs(error_bracket_1), std::fabs(error_s))};
  }
  double value_to_return =
      std::fabs(error_s) < std::fabs(error_bracket_1) ? s : a_bracket_1;
  otype_m->setGuess(&value_to_return);
}

template <class OptimizingClass>
double BrentsMethod<OptimizingClass>::calculateError(double* a_guess) {
  otype_m->setGuess(a_guess);
  return otype_m->calculateSignedScalarError();
}

}  // namespace IRL

#endif  // SRC_OPTIMIZATION_BRENTS_METHOD_TPP_
