// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_OPTIMIZATION_SECANT_H_
#define SRC_OPTIMIZATION_SECANT_H_

#include <cmath>
#include <utility>

#include "src/helpers/helper.h"

namespace IRL {
/// \brief A templated class driver that performs the Secant method
/// for optimization.
///
/// This class is a general Gauss Newton solver
/// that takes a pointer to a templated class
/// that must implement certain methods (detailed below).
/// Upon completion, the final solution should be stored
/// in the pointed-to-class because of the use of
/// `updateGuess()`.
///
/// Requirements for Optimizing class:
/// `calculateSignedScalarError()` : A method to calculate a scalar error
/// for the class with sign of (correct-guess).
/// - `updateGuess(...)` : A method that takes in the delta change
/// and computes a new guess vector (which it is storing itself)
/// - `updateBoundaries(...)` : (Only if `solveAndTrackBounds()` is being used)
/// A method to keep track of the best current solutions on the +/- side of the
/// zero.
/// - `errorTooHigh(...)` : A method that takes error
/// and returns a boolean whether the error is low
/// enough to stop optimization and return.
/// -`iterationTooHigh(...)` : A method that takes
/// the number of iterations and returns a bool
/// whether the maximum number of allowable iterations
/// has been exceeded.
template <class OptimizingClass>
class Secant {
 public:
  /// \brief Default constructor
  Secant(void) : otype_m(nullptr) {}

  /// \brief Set otype_m to a new pointer and then call `solve()`.
  ///
  /// \param[in] a_setup_otype Pointer to object that has
  /// been setup and will be optimized.
  /// \param[in] a_initial_delta Initial step to be taken.
  void solve(OptimizingClass* a_setup_otype, const double a_initial_delta);

  /// \brief Set pointer for otype_m and then call `solveAndTrackBounds()`
  void solveAndTrackBounds(OptimizingClass* a_setup_otype,
                           const double a_initial_delta);

  /// \brief Return reason for exiting by integer
  int getReason(void) const;

  /// \brief Default destructor
  ~Secant(void) = default;

 private:
  /// \brief Perform optimization.
  ///
  /// \param[in] a_initial_delta Initial step to be taken.
  void solve(const double a_initial_delta);

  /// \brief Perform optimization and track best answers on either side
  /// of the zero to finish by bisection if necessary.
  ///
  /// This requires the `updateBoundaries()` method mentioned in the
  /// class description.
  void solveAndTrackBounds(const double a_initial_delta);

  /// \brief Pointer to object of class `OptimizingClass`
  /// that is being optimized.
  OptimizingClass* otype_m;
  /// \brief Change in parameter being optimized.
  double delta_m;
  /// \brief Iterations of the Secant method.
  UnsignedIndex_t iteration_m;
  /// \brief Integer indicating reason for Secant method exiting.
  ///
  /// Reasons:
  /// - >= 0 : Number of iterations taken reduce error to acceptable level.
  /// - -1 : Exited due to exceeding maximum number of iterations.
  int reason_for_exit_m;
};

//******************************************************************* //
//     Templated function definitions placed below this
//******************************************************************* //

}  // namespace IRL

#include "src/optimization/secant.tpp"

#endif  // SRC_OPTIMIZATION_SECANT_H_
