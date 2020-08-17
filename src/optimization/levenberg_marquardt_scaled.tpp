// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_OPTIMIZATION_LEVENBERG_MARQUARDT_SCALED_TPP_
#define SRC_OPTIMIZATION_LEVENBERG_MARQUARDT_SCALED_TPP_

namespace IRL {

template <class OptimizingClass, int kColumns>
LevenbergMarquardtScaled<OptimizingClass, kColumns>::LevenbergMarquardtScaled(
    void)
    : otype_m(nullptr) {}

template <class OptimizingClass, int kColumns>
void LevenbergMarquardtScaled<OptimizingClass, kColumns>::solve(
    OptimizingClass *a_setup_otype, const int a_number_of_rows,
    const Eigen::Matrix<double, kColumns, 1> &a_jacobian_delta) {
  otype_m = a_setup_otype;
  this->solve(a_number_of_rows, a_jacobian_delta);
}

template <class OptimizingClass, int kColumns>
int LevenbergMarquardtScaled<OptimizingClass, kColumns>::getReason(void) {
  return reason_for_exit_m;
}

template <class OptimizingClass, int kColumns>
UnsignedIndex_t
LevenbergMarquardtScaled<OptimizingClass, kColumns>::getIterationCount(void) {
  return iteration_m;
}

template <class OptimizingClass, int kColumns>
void LevenbergMarquardtScaled<OptimizingClass, kColumns>::solve(
    const int a_number_of_rows,
    const Eigen::Matrix<double, kColumns, 1> &a_jacobian_delta) {
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
  Eigen::Matrix<double, Eigen::Dynamic, 1> reduction_step;
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
      // Exiting because minimum found (but error still abou acceptable_error)
      reason_for_exit_m = -2;
      return;
    }

    // Calculate new error and see if an improvement
    otype_m->updateGuess(&delta_m);
    double guess_error = otype_m->calculateScalarError();
    double actual_reduction = -1.0;
    if (0.1 * guess_error < error) {
      actual_reduction = 1.0 - std::pow(guess_error / error, 2.0);
    }

    // Predicted error reduction
    reduction_step = jacobian_transpose_m.transpose() * delta_m;

    const double temp1 = std::pow(reduction_step.stableNorm() / error, 2.0);
    // const double temp2 = lambda * std::pow(delta_m.stableNorm() /
    // error, 2.0);
    const double predicted_reduction = temp1; // + 2.0 * temp2;
    // const double directional_derivative = -(temp1 + temp2);

    // Ratio of actual to predicted
    const double ratio = predicted_reduction != 0.0
                             ? actual_reduction / predicted_reduction
                             : 0.0;

    // if (ratio <= 0.25) {
    //   double temp;
    //   if (actual_reduction >= 0.0) {
    //     temp = 0.5;
    //   } else {
    //     temp = 0.5 * directional_derivative /
    //            (directional_derivative + 0.5 * actual_reduction);
    //   }
    //   if (0.1 * guess_error >= error || temp < 0.1) {
    //     temp = 0.1;
    //   }
    //   lambda /= temp;
    // } else if (!(lambda != 0.0 && ratio < 0.75)) {
    //   lambda *= 0.5;
    // }
    // std::cout << "Reductions: " << actual_reduction << " " << temp1 << " "
    //          << temp2 << " " << predicted_reduction << " " << ratio
    //          << std::endl;
    // std::cout << "Reduction step " << reduction_step << std::endl;
    // std::cout << "Jacobian transpose " << jacobian_transpose_m <<
    // std::endl; std::cout << "Delta " << delta_m << std::endl;

    // If reconstruction is not an improvement, increase lambda and try again
    // Otherwise, accept solution and take another step
    if (ratio < 1.0e-4) {
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
void LevenbergMarquardtScaled<OptimizingClass, kColumns>::calculateJacobian(
    const Eigen::Matrix<double, kColumns, 1> &a_delta,
    Eigen::Matrix<double, kColumns, Eigen::Dynamic> *a_jacobian_transpose,
    Eigen::Matrix<double, kColumns, kColumns> *a_jacTjac) {
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

} // namespace IRL

#endif // SRC_OPTIMIZATION_LEVENBERG_MARQUARDT_SCALED_TPP_
