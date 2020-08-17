// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_OPTIMIZATION_LEVENBERG_MARQUARDT_H_
#define SRC_OPTIMIZATION_LEVENBERG_MARQUARDT_H_

#include <cmath>
#include <utility>

#include <Eigen/Dense>  // Eigen header

#include "src/helpers/helper.h"

namespace IRL {
/// \brief Levenberg-Marquardt optimization routine.
///
/// Requirements for OptimizingClass:
/// - `double calculateScalarError(void)` : A method to calculate a scalar error
/// that we are trying to minimize.
/// - `Eigen::Matrix<double,kRows,1>calculateVectorError(void)` : A method that
/// returns the vector return (correct_values - guess_values) by value
/// - `void updateGuess(Eigen::Matrix<double,kColumns,1>)` : A method that takes
/// in the delta change and computes a new guess vector (which it is storing
/// itself)
/// - `void updateBestGuess(void)` : A method that updates the best guess
/// and all other things necessary before a new Jacobian
/// is calculated and a new step is taken.
/// - `Eigen::Matrix<double,kRows,1>calculateChangeInGuess(void)` : A method
/// that calculates the difference between Guess variables and bestGuess (for
/// use in calculating derivative in Jacobian).
/// - void increaseLambda(double*) : A method to increase the value of
/// lambda (for failed attempts at finding a new minimum).
/// - void decreaseLambda(double*) : A method to decrease the value of
/// lambda (for successful attempts at finding a new minimum).
/// - `bool errorTooHigh(const double)` : A method that takes a scalar error
/// and returns a boolean whether the error is low
/// enough to stop optimization and return.
/// -`bool iterationTooHigh(const int)` : A method that takes
/// the number of iterations and returns a bool
/// whether the maximum number of allowable iterations
/// has been exceeded.
/// - `minimumReached(const Eigen::Matrix<double,kColumns,1>)` : A method that
/// takes delta and determines if the optimization has reached a minimum, return
/// a bool `true` if optimization should exit. `shouldComputeJacobian(const int,
/// const int)` : A method that returns a bool for whether or not a jacobian
/// should be computed when given the current iteration and the last iteration
/// the Jacobian was computed for.
///
/// kRows is the number of rows involved in the error vector of the
///  Levenberg-Marquardt system [y-f].
///
/// kColumns is the number of columns involved in the Jacobian for
/// the Levenberg-Marquardt system, equal to the number of parameters
/// being fit.
template <class OptimizingClass, int kRows, int kColumns>
class LevenbergMarquardt {
 public:
  /// \brief Default construction
  LevenbergMarquardt(void);

  /// \brief Assigns otype_m to a new pointer and then solves.
  void solve(OptimizingClass* a_setup_otype,
             const Eigen::Matrix<double, kColumns, 1>& a_jacobian_delta);

  /// \brief Return reason for exiting by integer
  int getReason(void);

  /// \brief Return the number of iterations it took until exit.
  UnsignedIndex_t getIterationCount(void);

  /// \brief Default dedstructor
  ~LevenbergMarquardt(void) = default;

 private:
  /// \brief Perform non-linear optimization.
  void solve(const Eigen::Matrix<double, kColumns, 1>& a_jacobian_delta);

  /// \brief Calculate jacobian using first-order finite difference.
  void calculateJacobian(
      const Eigen::Matrix<double, kColumns, 1>& a_delta,
      Eigen::Matrix<double, kColumns, kRows>* a_jacobian_tranpose,
      Eigen::Matrix<double, kColumns, kColumns>* a_jacTjac);

  /// \brief Pointer to object of class `OptimizingClass`
  /// that is being optimized.
  OptimizingClass* otype_m;
  /// \brief Iterations of Levenberg-Marquardt algorithm.
  UnsignedIndex_t iteration_m;
  /// \brief Integer indicating reason for Levenberg-Marquardt exiting.
  ///
  /// Reasons:
  /// - >= 0 : Number of iterations taken reduce error to acceptable level.
  /// - -1 : Exited due to exceeding maximum number of iterations.
  /// - -2 : Exited due to minimum reached (largest magnitude in delta
  /// less than set amount).
  int reason_for_exit_m;
  /// \brief Transpose of Jacobian for guess vector.
  Eigen::Matrix<double, kColumns, kRows> jacobian_transpose_m;
  /// \brief Matrix `jacTJac_m` =  `jacobian_transpose_m` *
  /// tranpose(jacobian_transpose_m).
  Eigen::Matrix<double, kColumns, kColumns> jacTjac_m;
  /// \brief Change in parameters being fit.
  Eigen::Matrix<double, kColumns, 1> delta_m;
  /// \brief (JacTJac_m + lambda*I), and preconditioned with Jacobi
  /// pre-conditioner.
  Eigen::Matrix<double, kColumns, kColumns> A_m;
  /// \brief Error vector of correct - guess
  Eigen::Matrix<double, kRows, 1> vector_error_m;
  /// \brief RHS of (JacTJac_m + lambda*I)*delta =
  /// `jacobian_transpose_m`*(`vector_error_m`)
  Eigen::Matrix<double, kColumns, 1> rhs_m;
  /// \brief `rhs_m` with Jacobi preconditioner applied.
  Eigen::Matrix<double, kColumns, 1> rhs_precond_m;
};

// Overload of LevenbergMarquardt class that requires run time setting of matrix
// rows.
template <class OptimizingClass, int kColumns>
class LevenbergMarquardt<OptimizingClass, -1, kColumns> {
 public:
  /// \brief Default construction
  LevenbergMarquardt(void);

  /// \brief Assigns otype_m to a new pointer and then solves.
  void solve(OptimizingClass* a_setup_otype, const int a_number_of_rows,
             const Eigen::Matrix<double, kColumns, 1>& a_jacobian_delta);

  /// \brief Return reason for exiting by integer
  int getReason(void);

  /// \brief Return the number of iterations it took until exit.
  UnsignedIndex_t getIterationCount(void);

  /// \brief Default dedstructor
  ~LevenbergMarquardt(void) = default;

 private:
  /// \brief Perform non-linear optimization to get
  /// construction according to R2P.
  void solve(const int a_number_of_rows,
             const Eigen::Matrix<double, kColumns, 1>& a_jacobian_delta);

  /// \brief Calculate jacobian using first-order finite difference.
  void calculateJacobian(
      const Eigen::Matrix<double, kColumns, 1>& a_delta,
      Eigen::Matrix<double, kColumns, Eigen::Dynamic>* a_jacobian_tranpose,
      Eigen::Matrix<double, kColumns, kColumns>* a_jacTjac);

  /// \brief Pointer to object of class `OptimizingClass`
  /// that is being optimized.
  OptimizingClass* otype_m;
  /// \brief Iterations of Levenberg-Marquardt algorithm.
  UnsignedIndex_t iteration_m;
  /// \brief Integer indicating reason for Levenberg-Marquardt exiting.
  ///
  /// Reasons:
  /// - >= 0 : Number of iterations taken reduce error to acceptable level.
  /// - -1 : Exited due to exceeding maximum number of iterations.
  /// - -2 : Exited due to minimum reached (largest magnitude in delta
  /// less than set amount).
  int reason_for_exit_m;
  /// \brief Transpose of Jacobian for guess vector.
  Eigen::Matrix<double, kColumns, Eigen::Dynamic> jacobian_transpose_m;
  /// \brief Matrix `jacTJac_m` =  `jacobian_transpose_m` *
  /// tranpose(jacobian_transpose_m).
  Eigen::Matrix<double, kColumns, kColumns> jacTjac_m;
  /// \brief Change in parameters being fit.
  Eigen::Matrix<double, kColumns, 1> delta_m;
  /// \brief (JacTJac_m + lambda*I), and preconditioned with Jacobi
  /// pre-conditioner.
  Eigen::Matrix<double, kColumns, kColumns> A_m;
  /// \brief Error vector of correct - guess
  Eigen::Matrix<double, Eigen::Dynamic, 1> vector_error_m;
  /// \brief RHS of (JacTJac_m + lambda*I)*delta =
  /// `jacobian_transpose_m`*(`vector_error_m`)
  Eigen::Matrix<double, kColumns, 1> rhs_m;
  /// \brief `rhs_m` with Jacobi preconditioner applied.
  Eigen::Matrix<double, kColumns, 1> rhs_precond_m;
};

}  // namespace IRL

#include "src/optimization/levenberg_marquardt.tpp"

#endif  // SRC_OPTIMIZATION_LEVENBERG_MARQUARDT_H_
