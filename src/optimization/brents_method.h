// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_OPTIMIZATION_BRENTS_METHOD_H_
#define SRC_OPTIMIZATION_BRENTS_METHOD_H_

#include <algorithm>
#include <cmath>
#include <utility>

#include "src/helpers/helper.h"

namespace IRL {
/// \brief A templated class driver that performs Brent's Method.
///
/// This class is a general Brent's Method solver
/// that takes a pointer to a templated class
/// that must implement certain methods (detailed below).
/// Upon completion, the final solution should be stored
/// in the pointed-to-class because of the use of
/// `setGuess()` before exiting `this->solve()`.
///
/// Requirements for Optimizing class:
/// `calculateSignedScalarError()` : A method to calculate a scalar error
/// for the class with sign of (correct-guess).
/// - `setGuess(...)` : A method that takes the guess value
/// and applies it to be used when `calculateSignedScalarError()` is called.
/// - `errorTooHigh(...)` : A method that takes error
/// and returns a boolean whether the error is low
/// enough to stop optimization and return.
template <class OptimizingClass>
class BrentsMethod {
 public:
  /// \brief Default constructor
  BrentsMethod(void);

  /// \brief Solution that resets optimizing class and brackets, then solves.
  ///
  /// \param[in] a_otype Pointer to class being optimized
  /// \param[in] a_bracket_0 Bracket on one side of zero.
  /// \param[in] a_bracket_1 Bracket on other side of zero from a_bracet_0.
  /// \param[in] a_delta Value for delta that dictates which method is used
  /// inside Brents Method.
  void solve(OptimizingClass* a_otype, double a_bracket_0, double a_bracket_1,
             const double a_delta);

  /// \brief Default destructor
  ~BrentsMethod(void) = default;

 private:
  /// \brief Given two locations that bracket a zero, find the zero.
  ///
  /// \param[in] a_bracket_0 Bracket on one side of zero.
  /// \param[in] a_bracket_1 Bracket on other side of zero from a_bracet_0.
  /// \param[in] a_delta Value for delta that dictates which method is used
  /// inside Brents Method.
  void solve(double a_bracket_0, double a_bracket_1, const double a_delta);

  /// \brief Set the guess and calculate the associated error.
  double calculateError(double* a_guess);

  /// \brief Pointer to object of class `OptimizingClass`
  /// that is being optimized.
  OptimizingClass* otype_m;
};

}  // namespace IRL

#include "src/optimization/brents_method.tpp"

#endif  // SRC_OPTIMIZATION_BRENTS_METHOD_H_
