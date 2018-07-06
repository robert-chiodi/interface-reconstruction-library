// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_OPTIMIZATION_LEVENBERG_MARQUARDT_TPP_
#define SRC_OPTIMIZATION_LEVENBERG_MARQUARDT_TPP_

namespace IRL {

template <class OptimizingClass, int kRows, int kColumns>
LevenbergMarquardt<OptimizingClass, kRows, kColumns>::LevenbergMarquardt(void)
    : otype_m(nullptr) {}

template <class OptimizingClass, int kRows, int kColumns>
void LevenbergMarquardt<OptimizingClass, kRows, kColumns>::solve(
    OptimizingClass* a_setup_otype,
    const Eigen::Matrix<double, kColumns, 1>& a_jacobian_delta) {
  otype_m = a_setup_otype;
  this->solve(a_jacobian_delta);
}

template <class OptimizingClass, int kRows, int kColumns>
int LevenbergMarquardt<OptimizingClass, kRows, kColumns>::getReason(void) {
  return reason_for_exit_m;
}

template <class OptimizingClass, int kRows, int kColumns>
UnsignedIndex_t
LevenbergMarquardt<OptimizingClass, kRows, kColumns>::getIterationCount(void) {
  return iteration_m;
}

template <class OptimizingClass, int kRows, int kColumns>
void LevenbergMarquardt<OptimizingClass, kRows, kColumns>::solve(
    const Eigen::Matrix<double, kColumns, 1>& a_jacobian_delta) {
  assert(otype_m != nullptr);
  // Calcualte initial error and save initial state
  delta_m = Eigen::Matrix<double, kColumns, 1>::Zero();
  otype_m->updateGuess(&delta_m);
  otype_m->updateBestGuess();

  double error = otype_m->calculateScalarError();
  reason_for_exit_m = 0;

  // Calculate initial Jacobian
  // This will change whatever guess variables the `OptimizingClass`
  // has, which is why calculateVectorError() is called before this.
  vector_error_m = otype_m->calculateVectorError();
  this->calculateJacobian(a_jacobian_delta, &jacobian_transpose_m, &jacTjac_m);

  // Calculate initial right hand side
  rhs_m = jacobian_transpose_m * vector_error_m;

  // Enter optimization loop
  double lambda = 1.0;
  iteration_m = 0;
  UnsignedIndex_t last_jacobian_iteration = 0;
  while (otype_m->errorTooHigh(error)) {
    iteration_m++;
    if (otype_m->iterationTooHigh(iteration_m)) {
      // Exiting because exceeding max iterations
      reason_for_exit_m = -1;
      return;
    }
    // Calculate A in A*delta = rhs
    A_m = jacTjac_m;
    for (int i = 0; i < kColumns; ++i) {
      A_m(i, i) += lambda;
      double jacobi_preconditioner = 1.0 / safelyEpsilon(A_m(i, i));
      for (int j = 0; j < kColumns; ++j) {
        A_m(i, j) *= jacobi_preconditioner;
      }
      rhs_precond_m(i) = rhs_m(i) * jacobi_preconditioner;
    }
    // Now solve system for delta
    delta_m = A_m.llt().solve(rhs_precond_m);

    // Check if delta is small, meaning minimum is reached
    if (otype_m->minimumReached(delta_m)) {
      // Exiting because minimum found (but error still about acceptable_error)
      reason_for_exit_m = -2;
      return;
    }

    // Calculate new error and see if an improvement
    otype_m->updateGuess(&delta_m);
    double guess_error = otype_m->calculateScalarError();

    // If reconstruction is not an improvement, increase lambda and try again
    // Otherwise, accept solution and take another step
    if (guess_error > error) {
      otype_m->increaseLambda(&lambda);
      continue;
    }

    // If accepted, decrease lambda to take a bigger step next time
    otype_m->decreaseLambda(&lambda);

    error = guess_error;
    otype_m->updateBestGuess();
    vector_error_m = otype_m->calculateVectorError();

    // Calculate Jacobian if enough iterations have past
    // This will change whatever guess variables the `OptimizingClass`
    // has, which is why calculateVectorError() is called before this.
    if (otype_m->shouldComputeJacobian(iteration_m, last_jacobian_iteration)) {
      this->calculateJacobian(a_jacobian_delta, &jacobian_transpose_m,
                              &jacTjac_m);
      last_jacobian_iteration = iteration_m;
    }

    // Get new right hand side using accepted (current best) reconstruction
    rhs_m = jacobian_transpose_m * vector_error_m;
  }
  // If exiting because error low enough, give number of iterations
  reason_for_exit_m = static_cast<int>(iteration_m);
}

template <class OptimizingClass, int kRows, int kColumns>
void LevenbergMarquardt<OptimizingClass, kRows, kColumns>::calculateJacobian(
    const Eigen::Matrix<double, kColumns, 1>& a_delta,
    Eigen::Matrix<double, kColumns, kRows>* a_jacobian_transpose,
    Eigen::Matrix<double, kColumns, kColumns>* a_jacTjac) {
  // Set up temporary delta
  Eigen::Matrix<double, kColumns, 1> solo_delta;
  // Calculate tranpose of Jacobian
  for (int parameter = 0; parameter < kColumns; ++parameter) {
    solo_delta = Eigen::Matrix<double, kColumns, 1>::Zero();
    solo_delta(parameter) = a_delta(parameter);
    otype_m->updateGuess(&solo_delta);
    Eigen::Matrix<double, kRows, 1> change_in_guess =
        otype_m->calculateChangeInGuess();
    for (int elem = 0; elem < kRows; ++elem) {
      (*a_jacobian_transpose)(parameter, elem) =
          change_in_guess(elem) / safelyEpsilon(solo_delta(parameter));
    }
  }
  // Calculate Tranpose(Jacobian) * Jacobian
  *a_jacTjac = (*a_jacobian_transpose) * (a_jacobian_transpose->transpose());
}

template <class OptimizingClass, int kColumns>
LevenbergMarquardt<OptimizingClass, -1, kColumns>::LevenbergMarquardt(void)
    : otype_m(nullptr) {}

template <class OptimizingClass, int kColumns>
void LevenbergMarquardt<OptimizingClass, -1, kColumns>::solve(
    OptimizingClass* a_setup_otype, const int a_number_of_rows,
    const Eigen::Matrix<double, kColumns, 1>& a_jacobian_delta) {
  otype_m = a_setup_otype;
  this->solve(a_number_of_rows, a_jacobian_delta);
}

template <class OptimizingClass, int kColumns>
int LevenbergMarquardt<OptimizingClass, -1, kColumns>::getReason(void) {
  return reason_for_exit_m;
}

template <class OptimizingClass, int kColumns>
UnsignedIndex_t
LevenbergMarquardt<OptimizingClass, -1, kColumns>::getIterationCount(void) {
  return iteration_m;
}

template <class OptimizingClass, int kColumns>
void LevenbergMarquardt<OptimizingClass, -1, kColumns>::solve(
    const int a_number_of_rows,
    const Eigen::Matrix<double, kColumns, 1>& a_jacobian_delta) {
  assert(otype_m != nullptr);

  // Construct actual matrices that are using Dynamic allocation
  jacobian_transpose_m = Eigen::Matrix<double, kColumns, Eigen::Dynamic>(
      kColumns, a_number_of_rows);
  vector_error_m =
      Eigen::Matrix<double, Eigen::Dynamic, 1>(a_number_of_rows, 1);

  // Calcualte initial error and save initial state
  delta_m = Eigen::Matrix<double, kColumns, 1>::Zero();
  otype_m->updateGuess(&delta_m);
  otype_m->updateBestGuess();

  double error = otype_m->calculateScalarError();
  reason_for_exit_m = 0;

  // Calculate initial Jacobian
  // This will change whatever guess variables the `OptimizingClass`
  // has, which is why calculateVectorError() is called before this.
  vector_error_m = otype_m->calculateVectorError();
  this->calculateJacobian(a_jacobian_delta, &jacobian_transpose_m, &jacTjac_m);

  // Calculate initial right hand side
  rhs_m = jacobian_transpose_m * vector_error_m;

  // Enter optimization loop
  double lambda = 1.0;
  iteration_m = 0;
  UnsignedIndex_t last_jacobian_iteration = 0;
  while (otype_m->errorTooHigh(error)) {
    iteration_m++;
    if (otype_m->iterationTooHigh(iteration_m)) {
      // Exiting because exceeding max iterations
      reason_for_exit_m = -1;
      return;
    }
    // Calculate A in A*delta = rhs
    A_m = jacTjac_m;
    for (int i = 0; i < kColumns; ++i) {
      A_m(i, i) += lambda;
      double jacobi_preconditioner = 1.0 / safelyEpsilon(A_m(i, i));
      for (int j = 0; j < kColumns; ++j) {
        A_m(i, j) *= jacobi_preconditioner;
      }
      rhs_precond_m(i) = rhs_m(i) * jacobi_preconditioner;
    }
    // Now solve system for delta
    delta_m = A_m.llt().solve(rhs_precond_m);

    // Check if delta is small, meaning minimum is reached
    if (otype_m->minimumReached(delta_m)) {
      // Exiting because minimum found (but error still about acceptable_error)
      reason_for_exit_m = -2;
      return;
    }

    // Calculate new error and see if an improvement
    otype_m->updateGuess(&delta_m);
    double guess_error = otype_m->calculateScalarError();

    // If reconstruction is not an improvement, increase lambda and try again
    // Otherwise, accept solution and take another step
    if (guess_error > error) {
      otype_m->increaseLambda(&lambda);
      continue;
    }

    // If accepted, decrease lambda to take a bigger step next time
    otype_m->decreaseLambda(&lambda);

    error = guess_error;
    otype_m->updateBestGuess();
    vector_error_m = otype_m->calculateVectorError();

    // Calculate Jacobian if enough iterations have past
    // This will change whatever guess variables the `OptimizingClass`
    // has, which is why calculateVectorError() is called before this.
    if (otype_m->shouldComputeJacobian(iteration_m, last_jacobian_iteration)) {
      this->calculateJacobian(a_jacobian_delta, &jacobian_transpose_m,
                              &jacTjac_m);
      last_jacobian_iteration = iteration_m;
    }

    // Get new right hand side using accepted (current best) reconstruction
    rhs_m = jacobian_transpose_m * vector_error_m;
  }
  // If exiting because error low enough, give number of iterations
  reason_for_exit_m = static_cast<int>(iteration_m);
}

template <class OptimizingClass, int kColumns>
void LevenbergMarquardt<OptimizingClass, -1, kColumns>::calculateJacobian(
    const Eigen::Matrix<double, kColumns, 1>& a_delta,
    Eigen::Matrix<double, kColumns, Eigen::Dynamic>* a_jacobian_transpose,
    Eigen::Matrix<double, kColumns, kColumns>* a_jacTjac) {
  // Set up temporary delta
  Eigen::Matrix<double, kColumns, 1> solo_delta;
  // Calculate tranpose of Jacobian
  for (int parameter = 0; parameter < kColumns; ++parameter) {
    solo_delta = Eigen::Matrix<double, kColumns, 1>::Zero();
    solo_delta(parameter) = a_delta(parameter);
    otype_m->updateGuess(&solo_delta);
    Eigen::Matrix<double, Eigen::Dynamic, 1> change_in_guess =
        otype_m->calculateChangeInGuess();
    for (int elem = 0; elem < change_in_guess.rows(); ++elem) {
      (*a_jacobian_transpose)(parameter, elem) =
          change_in_guess(elem) / safelyEpsilon(solo_delta(parameter));
    }
  }
  // Calculate Tranpose(Jacobian) * Jacobian
  *a_jacTjac = (*a_jacobian_transpose) * (a_jacobian_transpose->transpose());
}

}  // namespace IRL

#endif  // SRC_OPTIMIZATION_LEVENBERG_MARQUARDT_TPP_
