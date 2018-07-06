// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_OPTIMIZATION_BFGS_H_
#define SRC_OPTIMIZATION_BFGS_H_

#include <Eigen/Dense>  // Eigen header

#include "src/parameters/defined_types.h"
namespace IRL {

/// \brief Implementation of BFGS method according to the book
/// Numerical Optimization by Jorge Npcedal and Stephen Wright, 2006.
///
/// Working implementation of BFGS non-linear optimization algorithm. The line
/// search used is implemented below in the class
/// StrongWolfeConditionLineSearch. This line search preserves the Strong Wolfe
/// Conditions for step size and curvature, ensuring that the estimated Hessian
/// remains Symmetric Positive Definite and a valid search direction can always
/// be found.
///
/// Note: There is probably much opportunity for optimizing the implementation
/// of the algorithm. Right now it is rather slow.
template <class OptimizingClass, int kParameters>
class BFGS {
  friend class StrongWolfeConditionLineSearch;

  static constexpr double step_size_tolerance = 1.0e-10;

 public:
  BFGS(void);

  /// \brief Assigns otype_m to a new pointer and then solves.
  void solve(OptimizingClass* a_setup_otype,
             const Eigen::Matrix<double, kParameters, 1>& a_gradient_delta);

  /// \brief Return reason for exiting by integer
  int getReason(void);

  /// \brief Return the number of iterations it took until exit.
  UnsignedIndex_t getIterationCount(void);

  /// \brief Default dedstructor
  ~BFGS(void) = default;

 private:
  /// \brief Perform non-linear optimization.
  void solve(const Eigen::Matrix<double, kParameters, 1>& a_gradient_delta);

  /// \brief Calculate cost function with a given travel amount.
  /// This travel amount is the displacement (in parameter space)
  /// from the last accepted ('bestGuess') solution. If this is assumed to
  /// be represented as f(x) in typical BFGS literature, a travel amount
  /// would then be \alpha_k p_k, or the step distance times the
  /// step direction.
  double calculateCostFunction(
      Eigen::Matrix<double, kParameters, 1>* a_travel_amount);

  /// \brief Calculate cost function gradient.
  /// Calculated at a_starting_delta via forward difference with the step taken
  /// for finite-difference determined by parameter's entry in a_gradient_delta.
  Eigen::Matrix<double, kParameters, 1> calculateGradientOfCostFunction(
      const Eigen::Matrix<double, kParameters, 1>& a_gradient_delta);

  /// \brief Peform line search to determine step length
  Eigen::Matrix<double, kParameters, 1> getStep(
      const Eigen::Matrix<double, kParameters, 1>& a_search_direction);

  int reason_for_exit_m;
  UnsignedIndex_t iteration_m;
  OptimizingClass* otype_m;
};

/// \brief Implementation of line search algorithm from
/// "Line search algorithms with guaranteed sufficient decrease"
/// by Jorge More and David Thuente, ACM Transactions on Mathematical
/// Software, 1994. Some implementation details (such as the zoom function)
/// are also taken from the book Numerical Optimization by Jorge Npcedal and
/// Stephen Wright, 2006.
class StrongWolfeConditionLineSearch {
  static constexpr UnsignedIndex_t MAX_ITER = 10;
  // Constant to make sure we do not take too small a step size
  static constexpr double C1 = 1.0e-4;
  // Constant to make sure curvature of function is valid
  static constexpr double C2 = 0.9;
  // Step size to take when evaluating gradient along line
  static constexpr double GRADIENT_STEP = 1.0e-8;

 public:
  StrongWolfeConditionLineSearch(void) = default;

  // Requires optimizing class to implement a way to calculate \phi(\alpha) and
  // \phi'(\alpha), where \phi(\alpha) = f(x + \alpha p_k), with f being the
  // cost function we are searching for a minimum on and p_k being a search
  // direction (with p_k \cdot f' < 0). Returns the step length \alpha
  // alpha_max is the maximum allowable step size.
  template <class OptimizingClass, class EigenMatrix>
  EigenMatrix solve(OptimizingClass* a_otype,
                    const EigenMatrix& a_search_direction,
                    const double a_step_size_max);

 private:
  template <class OptimizingClass, class EigenMatrix>
  EigenMatrix zoom(OptimizingClass* a_otype,
                   const EigenMatrix& a_search_direction,
                   double a_step_size_lower, double a_step_size_upper,
                   double a_function_value_lower, double a_function_value_upper,
                   double a_gradient_lower, double a_gradient_upper);

  // Calculate the gradient of the cost function with respect to the step size.
  template <class OptimizingClass, class EigenMatrix>
  static double calculateGradientOfCostFunction(
      OptimizingClass* a_otype, const double a_step_size,
      const EigenMatrix& a_search_direction, const double a_step_size_delta);

  static double getNextTargetStep(
      const double a_step_size_lower, const double a_step_size_target,
      const double a_step_size_upper, const double a_function_value_lower,
      const double a_function_value_target, const double a_function_value_upper,
      const double a_gradient_lower, const double a_gradient_target,
      const double a_gradient_upper);

  static double cubicInterpolation(const double a_step_size_lower,
                                   const double a_step_size_upper,
                                   const double a_function_value_lower,
                                   const double a_function_value_upper,
                                   const double a_gradient_lower,
                                   const double a_gradient_upper);

  static double quadraticInterpolationWithTwoFunctions(
      const double a_step_size_lower, const double a_step_size_upper,
      const double a_function_value_lower, const double a_function_value_upper,
      const double a_gradient_lower);

  static double quadraticInterpolationWithTwoDerivatives(
      const double a_step_size_lower, const double a_step_size_upper,
      const double a_function_value_lower, const double a_gradient_lower,
      const double a_gradient_upper);

  double initial_function_value_m;
  double initial_function_gradient_m;
};

}  // namespace IRL

#include "src/optimization/bfgs.tpp"

#endif  // SRC_OPTIMIZATION_BFGS_H_
