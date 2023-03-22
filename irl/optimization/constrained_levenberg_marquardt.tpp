// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_OPTIMIZATION_CONSTRAINED_LEVENBERG_MARQUARDT_TPP_
#define IRL_OPTIMIZATION_CONSTRAINED_LEVENBERG_MARQUARDT_TPP_

namespace IRL {

template <class OptimizingClass, int kMeasures, int kParameters,
          int kConstraints>
ConstrainedLevenbergMarquardt<OptimizingClass, kMeasures, kParameters,
                              kConstraints>::ConstrainedLevenbergMarquardt(void)
    : otype_m(nullptr) {}

template <class OptimizingClass, int kMeasures, int kParameters,
          int kConstraints>
void ConstrainedLevenbergMarquardt<OptimizingClass, kMeasures, kParameters,
                                   kConstraints>::solve(OptimizingClass*
                                                            a_setup_otype) {
  otype_m = a_setup_otype;
  this->solve();
}

template <class OptimizingClass, int kMeasures, int kParameters,
          int kConstraints>
int ConstrainedLevenbergMarquardt<OptimizingClass, kMeasures, kParameters,
                                  kConstraints>::getReason(void) {
  return reason_for_exit_m;
}

template <class OptimizingClass, int kMeasures, int kParameters,
          int kConstraints>
UnsignedIndex_t
ConstrainedLevenbergMarquardt<OptimizingClass, kMeasures, kParameters,
                              kConstraints>::getIterationCount(void) {
  return iteration_m;
}

template <class OptimizingClass, int kMeasures, int kParameters,
          int kConstraints>
UnsignedIndex_t
ConstrainedLevenbergMarquardt<OptimizingClass, kMeasures, kParameters,
                              kConstraints>::getSubIterationCount(void) {
  return sub_iteration_m;
}

template <class OptimizingClass, int kMeasures, int kParameters,
          int kConstraints>
double ConstrainedLevenbergMarquardt<OptimizingClass, kMeasures, kParameters,
                                     kConstraints>::getConstraintError(void) {
  return std::sqrt(error_m);
}

template <class OptimizingClass, int kMeasures, int kParameters,
          int kConstraints>
Eigen::Matrix<double, kConstraints, 1>&
ConstrainedLevenbergMarquardt<OptimizingClass, kMeasures, kParameters,
                              kConstraints>::getConstraintErrorVector(void) {
  return constraints_m;
}

// Algorithm as given in http://www2.imm.dtu.dk/pubdb/edoc/imm5938.pdf with the
// augmented Lagrangian approach of
// http://www.seas.ucla.edu/~vandenbe/133B/lectures/nllseq.pdf
template <class OptimizingClass, int kMeasures, int kParameters,
          int kConstraints>
void ConstrainedLevenbergMarquardt<OptimizingClass, kMeasures, kParameters,
                                   kConstraints>::solve(void) {
  assert(otype_m != nullptr);

  // Initialize weights
  weights_m = otype_m->getWeights();

  // Initilize penalty parameter and Lagrange multipliers
  penalty_m = 1.0;
  multipliers_m = Eigen::Matrix<double, kConstraints, 1>::Zero();
  reason_for_exit_m = 0;

  Eigen::Matrix<double, kParameters, 1> final_delta;
  final_delta.setZero();

  // Calcualte initial error and save initial state
  delta_m = Eigen::Matrix<double, kParameters, 1>::Zero();
  otype_m->updateGuess(&delta_m, penalty_m, multipliers_m);
  otype_m->updateBestGuess();
  otype_m->storeBestGuess();
  otype_m->storeInitialGuess();

  // Initialize constraint errors
  error_m = DBL_MAX;
  double old_constraints_error = error_m;

  // Initialize constants for LM algorithm
  double initial_damping_factor = otype_m->getDampingParameterInitialFactor();
  int tries = 0;
  while (tries < 1 && otype_m->errorTooHigh(error_m)) {
    tries++;
    otype_m->restoreInitialGuess();

    // Enter optimization loop
    iteration_m = 0;
    sub_iteration_m = 0;

    // std::cout << "Starting augmented Lagrangian solver" << std::endl;
    while (otype_m->errorTooHigh(error_m)) {
      iteration_m++;
      if (otype_m->iterationTooHigh(iteration_m)) {
        // Exiting because exceeding max iterations
        reason_for_exit_m = -1;
        break;
      }

      // Calcualte initial error and save initial state
      delta_m.setZero();
      otype_m->updateGuess(&delta_m, penalty_m, multipliers_m);
      otype_m->updateBestGuess();

      // Calculate initial error vector and jacobian
      double error = otype_m->calculateScalarError(penalty_m, multipliers_m);
      otype_m->calculateVectorErrorAndJacobian(
          &vector_error_m, &jacobian_transpose_m, &constraints_m, penalty_m,
          multipliers_m);
      jacTjac_m =
          jacobian_transpose_m * weights_m * jacobian_transpose_m.transpose();

      // Calculate initial right hand side
      rhs_m = jacobian_transpose_m * weights_m * vector_error_m;

      // Enter sub-optimization loop
      double damping_param = std::numeric_limits<double>::lowest();
      double damping_factor = otype_m->getDampingParameterIncreaseFactor();
      double penalty_factor = otype_m->getPenaltyParameterIncreaseFactor();

      // The damping parameter is initialised proportional to the max diagonal
      // value of the jacobian...
      for (UnsignedIndex_t i = 0; i < kParameters; ++i) {
        damping_param = std::max(damping_param, jacTjac_m(i, i));
      }
      // ... then multiplied by the initial damping factor (typically a small
      // number)
      damping_param *= initial_damping_factor;
      // std::cout << "    damping_param = " << damping_param << std::endl;

      UnsignedIndex_t sub_iteration = 0;
      while (otype_m->maxErrorTooHigh(&rhs_m) &&
             !otype_m->subIterationTooHigh(sub_iteration)) {
        sub_iteration_m++;
        sub_iteration++;

        // Calculate A in A*delta = rhs
        A_m = jacTjac_m;
        damping_param = std::max(1.0e-7, damping_param);
        damping_param = std::min(1.0e7, damping_param);
        for (UnsignedIndex_t i = 0; i < kParameters; ++i) {
          A_m(i, i) += damping_param;
        }
        for (UnsignedIndex_t i = 0; i < kParameters; ++i) {
          // if (A_m(i, i) > 1.0e-6) {
          double precond = 1.0 / A_m(i, i);
          for (UnsignedIndex_t j = 0; j < kParameters; ++j) {
            A_m(i, j) *= precond;
          }
          rhs_normalized_m(i) = rhs_m(i) * penalty_m * precond;
          // } else {
          //   rhs_normalized_m(i) = rhs_m(i) * penalty_m;
          // }
        }
        // Now solve system for delta
        delta_m = A_m.completeOrthogonalDecomposition().solve(rhs_normalized_m);
        // Eigen::JacobiSVD<Eigen::MatrixXd> svd(A_m);
        // double cond = svd.singularValues()(0) /
        //               svd.singularValues()(svd.singularValues().size() - 1);
        // std::cout << "Matrix cond = " << cond << std::endl;

        // delta_m = (A_m.transpose() * A_m)
        //               .ldlt()
        //               .solve(A_m.transpose() * rhs_normalized_m);
        // double relative_error = (A_m * delta_m - rhs_normalized_m).norm() /
        //                         rhs_normalized_m.norm();  // norm() is L2
        //                         norm
        // std::cout << "The relative error is:\n" << relative_error <<
        // std::endl;
        otype_m->clipChange(&delta_m);
        // std::cout << "    delta = " << std::endl << delta_m << std::endl;

        // Check if delta is small, meaning minimum is reached
        if (!otype_m->subErrorTooHigh(delta_m.squaredNorm())) {
          break;
        }

        // Calculate new error and see if an improvement
        otype_m->updateGuess(&delta_m, penalty_m, multipliers_m);
        double guess_error =
            otype_m->calculateScalarError(penalty_m, multipliers_m);

        // Calculate rho = the "gain ratio"
        double gain_ratio = (error - guess_error) /
                            (delta_m.transpose() *
                             (damping_param * delta_m + rhs_m * penalty_m));

        // If gain ratio is negative, increase the damping parameter to get
        // closer to steepest descent and reduce step size, then return to start
        if (gain_ratio < 0.0) {
          damping_param *= damping_factor;
          damping_factor *= otype_m->getDampingParameterIncreaseFactor();
          error = guess_error;
          continue;
        }
        // else: we multiply by 2, if rho is small
        //       we divide by 3, if rho is large, so as to get closer to
        //       Gauss-Newton
        // The following implementation is a continuous blend of this:
        damping_param *=
            std::max(1.0 / 3.0, (1.0 - std::pow(2.0 * gain_ratio - 1.0, 3.0)));
        damping_factor = otype_m->getDampingParameterIncreaseFactor();

        // // Backtracking line search
        // double alpha = 2.0;
        // double c = 0.25;
        // double tau = 0.5;
        // double t = c * 2.0 * rhs_m.transpose() * penalty_m * delta_m;
        // Eigen::Matrix<double, kParameters, 1> new_delta;
        // new_delta = alpha * delta_m;
        // otype_m->updateGuess(&new_delta, penalty_m, multipliers_m);
        // guess_error = otype_m->calculateScalarError(penalty_m,
        // multipliers_m); int iter_armijo = 0; while (guess_error > error +
        // alpha * t && iter_armijo < 50) {
        //   alpha *= tau;
        //   new_delta = alpha * delta_m;
        //   otype_m->updateGuess(&new_delta, penalty_m, multipliers_m);
        //   guess_error = otype_m->calculateScalarError(penalty_m,
        //   multipliers_m); iter_armijo++;
        // }
        // if (iter_armijo > 0)
        //   std::cout << " Armijo rule iterations = " << iter_armijo <<
        //   std::endl;
        // delta_m = new_delta;
        error = guess_error;

        final_delta += delta_m;
        // Store the current guess and update errors + jacobian
        otype_m->updateBestGuess();
        otype_m->calculateVectorErrorAndJacobian(
            &vector_error_m, &jacobian_transpose_m, &constraints_m, penalty_m,
            multipliers_m);
        // Compute new constraint error
        // error_m = 0.0;
        // for (UnsignedIndex_t i = 0; i < kConstraints; ++i) {
        //   const double error_i =
        //       std::min(-constraints_m(i), multipliers_m(i) / (2.0 *
        //       penalty_m));
        //   error_m += error_i * error_i;
        // }
        error_m = constraints_m.squaredNorm();
        // if (error_m < old_constraints_error) {
        //   otype_m->storeBestGuess();
        // }
        jacTjac_m =
            jacobian_transpose_m * weights_m * jacobian_transpose_m.transpose();

        // Get new right hand side using accepted (current best) reconstruction
        rhs_m = jacobian_transpose_m * weights_m * vector_error_m / penalty_m;
      }

      // if (error_m > old_constraints_error || error_m > 1.0e-2) {
      //   otype_m->restoreBestGuess();
      //   continue;
      // }

      // Update multipliers and penalty parameter
      if (iteration_m > 5 || !otype_m->subIterationTooHigh(sub_iteration)) {
        multipliers_m += 2.0 * penalty_m * constraints_m;
        // error_m = std::max(error_m, constraints_m.squaredNorm());
        if (iteration_m > 1 && error_m >= 0.25 * 0.25 * old_constraints_error) {
          penalty_m *= penalty_factor;
        }
      }

      for (UnsignedIndex_t i = 0; i < kConstraints; ++i) {
        multipliers_m(i) = std::max(multipliers_m(i), -1.0e4);
        multipliers_m(i) = std::min(multipliers_m(i), 1.0e4);
      }

      // Store constraint error
      old_constraints_error = error_m;
      error_m = constraints_m.squaredNorm();

      // std::cout << "  Iter : " << iteration_m
      //           << "  Sub-Iter : " << sub_iteration_m
      //           << " -- Penalty = " << penalty_m
      //           << " -- Multiplier 0 = " << multipliers_m(0)
      //           << " -- Multiplier 1 = " << multipliers_m(1)
      //           << " -- Multiplier 2 = " << multipliers_m(2)
      //           << " -- Multiplier 3 = " << multipliers_m(3)
      //           << " -- Error = " << std::sqrt(error_m) << std::endl;

      // std::cout << "   Final delta = " << std::endl << final_delta <<
      // std::endl; If exiting because error low enough, give number of
      // iterations
      reason_for_exit_m = static_cast<int>(iteration_m);

      if (initial_damping_factor > 1.0e-6) initial_damping_factor /= 5.0;
    }
    // otype_m->restoreBestGuess();
  }
}
}  // namespace IRL

#endif  // IRL_OPTIMIZATION_CONSTRAINED_LEVENBERG_MARQUARDT_TPP_
