// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_OPTIMIZATION_BFGS_TPP_
#define SRC_OPTIMIZATION_BFGS_TPP_

#include <float.h>
#include <array>
#include <cassert>
#include <cmath>

#include <iostream>

#include "src/helpers/helper.h"

namespace IRL {

template <class OptimizingClass, int kParameters>
BFGS<OptimizingClass, kParameters>::BFGS(void) : otype_m{nullptr} {}

template <class OptimizingClass, int kParameters>
void BFGS<OptimizingClass, kParameters>::solve(
    OptimizingClass* a_setup_otype,
    const Eigen::Matrix<double, kParameters, 1>& a_gradient_delta) {
  otype_m = a_setup_otype;
  this->solve(a_gradient_delta);
}

template <class OptimizingClass, int kParameters>
int BFGS<OptimizingClass, kParameters>::getReason(void) {
  return reason_for_exit_m;
}

template <class OptimizingClass, int kParameters>
UnsignedIndex_t BFGS<OptimizingClass, kParameters>::getIterationCount(void) {
  return iteration_m;
}

template <class OptimizingClass, int kParameters>
void BFGS<OptimizingClass, kParameters>::solve(
    const Eigen::Matrix<double, kParameters, 1>& a_gradient_delta) {
  assert(otype_m != nullptr);

  // Calculate initial error and save initial state
  Eigen::Matrix<double, kParameters, 1> step_in_params =
      Eigen::Matrix<double, kParameters, 1>::Zero();
  otype_m->updateGuess(&step_in_params);
  otype_m->updateBestGuess();
  double error = otype_m->calculateScalarError();
  Eigen::Matrix<double, kParameters, 1> gradient_of_cost_function =
      this->calculateGradientOfCostFunction(a_gradient_delta);
  Eigen::Matrix<double, kParameters, 1> old_gradient_of_cost_function =
      gradient_of_cost_function;

  // Set initial inverse Hessian matrix as identity
  Eigen::Matrix<double, kParameters, kParameters> identity_matrix =
      Eigen::Matrix<double, kParameters, kParameters>::Identity(kParameters,
                                                                kParameters);
  Eigen::Matrix<double, kParameters, kParameters> inv_Hessian = identity_matrix;

  iteration_m = 0;
  reason_for_exit_m = 0;
  Eigen::Matrix<double, kParameters, 1> search_direction;
  Eigen::Matrix<double, kParameters, 1> change_in_gradient;

  while (otype_m->errorTooHigh(error)) {
    // Compute search direction
    search_direction = -inv_Hessian * gradient_of_cost_function;

    // Find step through a line search. Make sure step satisfies strong Wolfe
    // conditions

    // Calculate step in parameters and change in gradient
    // Step is calculated with a line search that obeys Strong Wolfe Conditions
    // What is returned is the vector of displacements for each parameter.
    // This is done because the step might be slightly perturbed from the given
    // search direction since some OptimizingClasses may impose constraints,
    // slightly altering the displacements given to them (such as R2P
    // modifying change in distances in order to maintain volume conservation).
    step_in_params = this->getStep(search_direction);

    // Do not need to re-update guess. Last performed update should
    // have been to get this step_in_params.
    // otype_m->updateGuess(&step_in_params);
    otype_m->updateBestGuess();
    error = otype_m->calculateScalarError();

    gradient_of_cost_function =
        this->calculateGradientOfCostFunction(a_gradient_delta);

    change_in_gradient =
        gradient_of_cost_function - old_gradient_of_cost_function;

    // Update size of H if first iteration to prevent erroneous first large step
    // from the initial guess of the Identity matrix for the inverse Hessian.
    if (iteration_m == 0) {
      double scaling = (change_in_gradient.transpose() * step_in_params);
      scaling /= (change_in_gradient.transpose() * change_in_gradient);
      inv_Hessian = scaling * identity_matrix;
    }

    // Update approximation of inverse Hessian
    const double rho_factor =
        1.0 / (change_in_gradient.transpose() * step_in_params);
    inv_Hessian = (identity_matrix - rho_factor * step_in_params *
                                         change_in_gradient.transpose()) *
                      inv_Hessian *
                      (identity_matrix - rho_factor * change_in_gradient *
                                             step_in_params.transpose()) +
                  rho_factor * step_in_params * step_in_params.transpose();

    // Check if minimum reached
    if (otype_m->minimumReached(gradient_of_cost_function)) {
      // Exiting because minimum found (but error still about acceptable_error)
      reason_for_exit_m = -2;
      return;
    }

    // Check if minimum reached with very small step.
    if (step_in_params.squaredNorm() < BFGS::step_size_tolerance) {
      reason_for_exit_m = -2;
      return;
    }

    // Check if max iteration exceeded.
    if (otype_m->iterationTooHigh(iteration_m)) {
      reason_for_exit_m = -1;
      return;
    }

    old_gradient_of_cost_function = gradient_of_cost_function;
    iteration_m++;
  }

  // If exiting because error low enough, give number of iterations
  reason_for_exit_m = static_cast<int>(iteration_m);
}

template <class OptimizingClass, int kParameters>
Eigen::Matrix<double, kParameters, 1>
BFGS<OptimizingClass, kParameters>::getStep(
    const Eigen::Matrix<double, kParameters, 1>& a_search_direction) {
  StrongWolfeConditionLineSearch line_search;
  static constexpr double max_step_size = 1.5;
  return line_search.solve(this, a_search_direction, max_step_size);
}

template <class OptimizingClass, int kParameters>
double BFGS<OptimizingClass, kParameters>::calculateCostFunction(
    Eigen::Matrix<double, kParameters, 1>* a_travel_amount) {
  otype_m->updateGuess(a_travel_amount);
  return otype_m->calculateScalarError();
}

template <class OptimizingClass, int kParameters>
Eigen::Matrix<double, kParameters, 1>
BFGS<OptimizingClass, kParameters>::calculateGradientOfCostFunction(
    const Eigen::Matrix<double, kParameters, 1>& a_gradient_delta) {
  double starting_error = otype_m->calculateScalarError();
  Eigen::Matrix<double, kParameters, 1> gradient_of_cost_function;
  Eigen::Matrix<double, kParameters, 1> solo_delta;
  for (int parameter = 0; parameter < kParameters; ++parameter) {
    solo_delta = Eigen::Matrix<double, kParameters, 1>::Zero();
    solo_delta(parameter) = a_gradient_delta(parameter);
    const double current_error = this->calculateCostFunction(&solo_delta);
    gradient_of_cost_function(parameter) =
        (current_error - starting_error) / safelyEpsilon(solo_delta(parameter));
  }
  return gradient_of_cost_function;
}

template <class OptimizingClass, class EigenMatrix>
EigenMatrix StrongWolfeConditionLineSearch::solve(
    OptimizingClass* a_otype, const EigenMatrix& a_search_direction,
    const double a_step_size_max) {
  assert(a_otype != nullptr);
  assert(a_step_size_max > 0.0);

  EigenMatrix parameter_displacement;

  // Select initial step size
  double old_step_size = 0.0;
  double step_size = a_step_size_max > 1.0 ? 1.0 : 0.5 * a_step_size_max;

  // Calculate beginning cost function value and gradient.
  parameter_displacement = EigenMatrix::Zero();
  initial_function_value_m =
      a_otype->calculateCostFunction(&parameter_displacement);
  initial_function_gradient_m =
      StrongWolfeConditionLineSearch::calculateGradientOfCostFunction(
          a_otype, 0.0, a_search_direction, GRADIENT_STEP);
  double old_function_value = initial_function_value_m;
  double old_function_gradient = initial_function_gradient_m;

  UnsignedIndex_t iteration = 1;
  while (iteration < StrongWolfeConditionLineSearch::MAX_ITER) {
    parameter_displacement = step_size * a_search_direction;
    double current_function_value =
        a_otype->calculateCostFunction(&parameter_displacement);
    double current_function_gradient =
        StrongWolfeConditionLineSearch::calculateGradientOfCostFunction(
            a_otype, step_size, a_search_direction, GRADIENT_STEP);

    if ((current_function_value >
         initial_function_value_m + StrongWolfeConditionLineSearch::C1 *
                                        step_size *
                                        initial_function_gradient_m) ||
        (current_function_value >= old_function_value && iteration > 1)) {
      return StrongWolfeConditionLineSearch::zoom(
          a_otype, a_search_direction, old_step_size, step_size,
          old_function_value, current_function_value, old_function_gradient,
          current_function_gradient);
    }
    if (fabs(current_function_gradient) <=
        -StrongWolfeConditionLineSearch::C2 * initial_function_gradient_m) {
      return parameter_displacement;
    }
    if (current_function_gradient >= 0.0) {
      return StrongWolfeConditionLineSearch::zoom(
          a_otype, a_search_direction, step_size, old_step_size,
          current_function_value, old_function_value, current_function_gradient,
          old_function_gradient);
    }

    old_step_size = step_size;
    old_function_value = current_function_value;
    old_function_gradient = current_function_gradient;
    step_size *= 1.2;
    ++iteration;
  }
  return parameter_displacement;  // Hope this is good enough if exceeded max
                                  // iterations..
}
template <class OptimizingClass, class EigenMatrix>
EigenMatrix StrongWolfeConditionLineSearch::zoom(
    OptimizingClass* a_otype, const EigenMatrix& a_search_direction,
    double a_step_size_lower, double a_step_size_upper,
    double a_function_value_lower, double a_function_value_upper,
    double a_gradient_lower, double a_gradient_upper) {
  double step_size = 0.5 * (a_step_size_upper + a_step_size_lower);
  std::array<double, 3> bracket_width{a_step_size_upper - a_step_size_lower,
                                      DBL_MAX, DBL_MAX};
  EigenMatrix parameter_displacement;
  UnsignedIndex_t iteration = 0;
  while (iteration < StrongWolfeConditionLineSearch::MAX_ITER) {
    parameter_displacement = step_size * a_search_direction;
    double current_function_value =
        a_otype->calculateCostFunction(&parameter_displacement);
    double current_function_gradient =
        StrongWolfeConditionLineSearch::calculateGradientOfCostFunction(
            a_otype, step_size, a_search_direction, GRADIENT_STEP);
    double new_step_size = StrongWolfeConditionLineSearch::getNextTargetStep(
        a_step_size_lower, step_size, a_step_size_upper, a_function_value_lower,
        current_function_value, a_function_value_upper, a_gradient_lower,
        current_function_gradient, a_gradient_upper);

    if (current_function_value >
            initial_function_value_m + StrongWolfeConditionLineSearch::C1 *
                                           step_size *
                                           initial_function_gradient_m ||
        current_function_value >= a_function_value_lower) {
      a_step_size_upper = step_size;
      a_function_value_upper = current_function_value;
      a_gradient_upper = current_function_gradient;
    } else {
      if (fabs(current_function_gradient) <=
          -StrongWolfeConditionLineSearch::C2 * initial_function_gradient_m) {
        return parameter_displacement;
      }
      if (current_function_gradient * (a_step_size_upper - a_step_size_lower) >=
          0.0) {
        a_step_size_upper = a_step_size_lower;
        a_function_value_upper = a_function_value_lower;
        a_gradient_upper = a_gradient_lower;
      }
      a_step_size_lower = step_size;
      a_function_value_lower = current_function_value;
      a_gradient_lower = current_function_gradient;
    }

    step_size = new_step_size;
    bracket_width[2] = bracket_width[1];
    bracket_width[1] = bracket_width[0];
    bracket_width[0] = a_step_size_upper - a_step_size_lower;
    if (bracket_width[0] < 10.0 * DBL_EPSILON) {
      break;
    }
    static constexpr double delta_for_bisect = 2.0 / 3.0;
    if (iteration > 1 &&
        (bracket_width[0] / bracket_width[2]) > delta_for_bisect) {
      step_size = a_step_size_lower + 0.5 * bracket_width[0];
    }

    ++iteration;
  }
  return parameter_displacement;
}

template <class OptimizingClass, class EigenMatrix>
double StrongWolfeConditionLineSearch::calculateGradientOfCostFunction(
    OptimizingClass* a_otype, const double a_step_size,
    const EigenMatrix& a_search_direction, const double a_step_size_delta) {
  EigenMatrix gradient_location = a_step_size * a_search_direction;
  EigenMatrix perturbed_location =
      (a_step_size + a_step_size_delta) * a_search_direction;
  double step_size_function_value =
      a_otype->calculateCostFunction(&gradient_location);
  double delta_moved_function_value =
      a_otype->calculateCostFunction(&perturbed_location);
  return (delta_moved_function_value - step_size_function_value) /
         safelyTiny((perturbed_location - gradient_location).norm());
}

inline double StrongWolfeConditionLineSearch::getNextTargetStep(
    const double a_step_size_lower, const double a_step_size_target,
    const double a_step_size_upper, const double a_function_value_lower,
    const double a_function_value_target, const double a_function_value_upper,
    const double a_gradient_lower, const double a_gradient_target,
    const double a_gradient_upper) {
  double step_to_return = -DBL_MAX;  // Should never stay this vaulue;
  if (a_function_value_target > a_function_value_lower) {
    const double cubic_interpolant =
        StrongWolfeConditionLineSearch::cubicInterpolation(
            a_step_size_lower, a_step_size_target, a_function_value_lower,
            a_function_value_target, a_gradient_lower, a_gradient_target);

    const double quadratic_interpolant =
        StrongWolfeConditionLineSearch::quadraticInterpolationWithTwoFunctions(
            a_step_size_lower, a_step_size_target, a_function_value_lower,
            a_function_value_target, a_gradient_lower);

    step_to_return = fabs(cubic_interpolant - a_step_size_lower) <
                             fabs(quadratic_interpolant - a_step_size_lower)
                         ? cubic_interpolant
                         : 0.5 * (quadratic_interpolant + cubic_interpolant);

  } else if (a_function_value_target <= a_function_value_lower &&
             a_gradient_lower * a_gradient_target < 0.0) {
    const double cubic_interpolant =
        StrongWolfeConditionLineSearch::cubicInterpolation(
            a_step_size_lower, a_step_size_target, a_function_value_lower,
            a_function_value_target, a_gradient_lower, a_gradient_target);

    const double quadratic_interpolant = StrongWolfeConditionLineSearch::
        quadraticInterpolationWithTwoDerivatives(
            a_step_size_lower, a_step_size_target, a_function_value_lower,
            a_gradient_lower, a_gradient_target);

    step_to_return = fabs(cubic_interpolant - a_step_size_target) >=
                             fabs(quadratic_interpolant - a_step_size_target)
                         ? cubic_interpolant
                         : quadratic_interpolant;

  } else if (a_function_value_target <= a_function_value_lower &&
             a_gradient_lower * a_gradient_target >= 0.0 &&
             fabs(a_gradient_target) <= fabs(a_gradient_lower)) {
    const double cubic_interpolant =
        StrongWolfeConditionLineSearch::cubicInterpolation(
            a_step_size_lower, a_step_size_target, a_function_value_lower,
            a_function_value_target, a_gradient_lower, a_gradient_target);

    const double quadratic_interpolant = StrongWolfeConditionLineSearch::
        quadraticInterpolationWithTwoDerivatives(
            a_step_size_lower, a_step_size_target, a_function_value_lower,
            a_gradient_lower, a_gradient_target);

    step_to_return = fabs(cubic_interpolant - a_step_size_target) <
                             fabs(quadratic_interpolant - a_step_size_target)
                         ? cubic_interpolant
                         : quadratic_interpolant;
    static constexpr double damping_delta_from_paper = 2.0 / 3.0;
    const double fall_back_step =
        a_step_size_target +
        damping_delta_from_paper * (a_step_size_upper - a_step_size_target);
    step_to_return = a_step_size_target > a_step_size_lower
                         ? std::min(fall_back_step, step_to_return)
                         : std::max(fall_back_step, step_to_return);

  } else if (a_function_value_target <= a_function_value_lower &&
             a_gradient_lower * a_gradient_target >= 0.0 &&
             fabs(a_gradient_target) > fabs(a_gradient_lower)) {
    const double cubic_interpolant =
        StrongWolfeConditionLineSearch::cubicInterpolation(
            a_step_size_upper, a_step_size_target, a_function_value_upper,
            a_function_value_target, a_gradient_upper, a_gradient_target);

    step_to_return = cubic_interpolant;

  } else {
    std::cout << "Should not be able to reach this. Debugging statement"
              << std::endl;
    std::exit(-1);
    return -DBL_MAX;  // should never happen
  }

  return step_to_return;
}

inline double StrongWolfeConditionLineSearch::cubicInterpolation(
    const double a_step_size_lower, const double a_step_size_upper,
    const double a_function_value_lower, const double a_function_value_upper,
    const double a_gradient_lower, const double a_gradient_upper) {
  const double d1 = a_gradient_lower + a_gradient_upper -
                    3.0 * (a_function_value_lower - a_function_value_upper) /
                        (a_step_size_lower - a_step_size_upper);
  const double d2 = std::copysign(
      std::sqrt(fabs(d1 * d1 - a_gradient_lower * a_gradient_upper)),
      a_step_size_upper - a_step_size_lower);

  double new_step =
      a_step_size_upper - (a_step_size_upper - a_step_size_lower) *
                              (a_gradient_upper + d2 - d1) /
                              (a_gradient_upper - a_gradient_lower + 2.0 * d2);

  return new_step;
}

inline double
StrongWolfeConditionLineSearch::quadraticInterpolationWithTwoFunctions(
    const double a_step_size_lower, const double a_step_size_upper,
    const double a_function_value_lower, const double a_function_value_upper,
    const double a_gradient_lower) {
  return a_step_size_lower -
         0.5 * ((a_step_size_lower - a_step_size_upper) * a_gradient_lower) /
             (a_gradient_lower -
              (a_function_value_lower - a_function_value_upper) /
                  (a_step_size_lower - a_step_size_upper));
}

inline double
StrongWolfeConditionLineSearch::quadraticInterpolationWithTwoDerivatives(
    const double a_step_size_lower, const double a_step_size_upper,
    const double a_function_value_lower, const double a_gradient_lower,
    const double a_gradient_upper) {
  return a_step_size_lower - (a_step_size_lower - a_step_size_upper) /
                                 (a_gradient_lower - a_gradient_upper) *
                                 a_gradient_lower;
}

}  // namespace IRL

#endif  // SRC_OPTIMIZATION_BFGS_TPP_
