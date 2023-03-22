// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_OPTIMIZATION_CONSTRAINED_LEVENBERG_MARQUARDT_H_
#define IRL_OPTIMIZATION_CONSTRAINED_LEVENBERG_MARQUARDT_H_

#include <cmath>
#include <utility>

#include <Eigen/Dense>  // Eigen header

#include "irl/helpers/helper.h"

namespace IRL {
/// \brief Levenberg-Marquardt optimization routine for constrained non-linear
/// least-squares problems
template <class OptimizingClass, int kMeasures, int kParameters,
          int kConstraints>
class ConstrainedLevenbergMarquardt {
 public:
  /// \brief Default construction
  ConstrainedLevenbergMarquardt(void);

  /// \brief Assigns otype_m to a new pointer and then solves using exact
  /// gradients.
  void solve(OptimizingClass* a_setup_otype);

  /// \brief Return reason for exiting by integer
  int getReason(void);

  /// \brief Return the number of iterations it took until exit.
  UnsignedIndex_t getIterationCount(void);

  /// \brief Return the number of sub-iterations it took until exit.
  UnsignedIndex_t getSubIterationCount(void);

  /// \brief Return the number of sub-iterations it took until exit.
  double getConstraintError(void);
  Eigen::Matrix<double, kConstraints, 1>& getConstraintErrorVector(void);

  /// \brief Default dedstructor
  ~ConstrainedLevenbergMarquardt(void) = default;

 private:
  /// \brief Perform non-linear optimization using exact gradients.
  void solve(void);

  /// \brief Pointer to object of class `OptimizingClass`
  /// that is being optimized.
  OptimizingClass* otype_m;
  /// \brief Iterations of augmented Lagrangian algorithm.
  UnsignedIndex_t iteration_m;
  /// \brief Iterations of Levenberg-Marquardt algorithm.
  UnsignedIndex_t sub_iteration_m;
  /// \brief Integer indicating reason for Levenberg-Marquardt exiting.
  ///
  /// Reasons:
  /// - >= 0 : Number of iterations taken reduce error to acceptable level.
  /// - -1 : Exited due to exceeding maximum number of iterations.
  int reason_for_exit_m;
  /// \brief Transpose of Jacobian for guess vector.
  Eigen::Matrix<double, kParameters, kMeasures + kConstraints>
      jacobian_transpose_m;
  /// \brief Matrix `jacTJac_m` =  `jacobian_transpose_m` *
  /// tranpose(jacobian_transpose_m).
  Eigen::Matrix<double, kParameters, kParameters> jacTjac_m;
  /// \brief Change in parameters being fit.
  Eigen::Matrix<double, kParameters, 1> delta_m;
  /// \brief (JacTJac_m + lambda*I), and preconditioned with Jacobi
  /// pre-conditioner.
  Eigen::Matrix<double, kParameters, kParameters> A_m;
  /// \brief Error vector of correct - guess
  Eigen::Matrix<double, kMeasures + kConstraints, 1> vector_error_m;
  /// \brief Constraints error vector.
  Eigen::Matrix<double, kConstraints, 1> constraints_m;
  /// \brief Array of weights.
  Eigen::DiagonalMatrix<double, kMeasures + kConstraints> weights_m;
  /// \brief RHS of (JacTJac_m + lambda*I)*delta =
  /// `jacobian_transpose_m`*(`vector_error_m`)
  Eigen::Matrix<double, kParameters, 1> rhs_m;
  Eigen::Matrix<double, kParameters, 1> rhs_normalized_m;
  /// \brief Array of Lagrange multipliers.
  Eigen::Matrix<double, kConstraints, 1> multipliers_m;
  /// \brief Penalty parameter.
  double penalty_m;
  /// \brief Final error.
  double error_m;
};

// // Overload of LevenbergMarquardt class that requires run time setting of
// matrix
// // rows.
// template <class OptimizingClass, int kMeasures, int kParameters,
//           int kConstraints>
// class ConstrainedLevenbergMarquardt<OptimizingClass, -1, kParameters,
//                                     kConstraints> {
//  public:
//   /// \brief Default construction
//   ConstrainedLevenbergMarquardt(void);

//   /// \brief Assigns otype_m to a new pointer and then solves using exact
//   /// gradients.
//   void solve(OptimizingClass* a_setup_otype);

//   /// \brief Return reason for exiting by integer
//   int getReason(void);

//   /// \brief Return the number of iterations it took until exit.
//   UnsignedIndex_t getIterationCount(void);

//   /// \brief Default dedstructor
//   ~ConstrainedLevenbergMarquardt(void) = default;

//  private:
//   /// \brief Perform non-linear optimization using exact gradients.
//   void solve(void);

//   /// \brief Pointer to object of class `OptimizingClass`
//   /// that is being optimized.
//   OptimizingClass* otype_m;
//   /// \brief Iterations of Levenberg-Marquardt algorithm.
//   UnsignedIndex_t iteration_m;
//   /// \brief Integer indicating reason for Levenberg-Marquardt exiting.
//   ///
//   /// Reasons:
//   /// - >= 0 : Number of iterations taken reduce error to acceptable level.
//   /// - -1 : Exited due to exceeding maximum number of iterations.
//   /// - -2 : Exited due to minimum reached (largest magnitude in delta
//   /// less than set amount).
//   int reason_for_exit_m;
//   /// \brief Transpose of Jacobian for guess vector.
//   Eigen::Matrix<double, kParameters + kConstraints, Eigen::Dynamic>
//       jacobian_transpose_m;
//   /// \brief Matrix `jacTJac_m` =  `jacobian_transpose_m` *
//   /// tranpose(jacobian_transpose_m).
//   Eigen::Matrix<double, kParameters + kConstraints, kParameters +
//   kConstraints>
//       jacTjac_m;
//   /// \brief Change in parameters being fit.
//   Eigen::Matrix<double, kParameters + kConstraints, 1> delta_m;
//   /// \brief (JacTJac_m + lambda*I), and preconditioned with Jacobi
//   /// pre-conditioner.
//   Eigen::Matrix<double, kParameters + kConstraints, kParameters +
//   kConstraints>
//       A_m;
//   /// \brief Error vector of correct - guess
//   Eigen::Matrix<double, Eigen::Dynamic, 1> vector_error_m;
//   /// \brief RHS of (JacTJac_m + lambda*I)*delta =
//   /// `jacobian_transpose_m`*(`vector_error_m`)
//   Eigen::Matrix<double, kParameters + kConstraints, 1> rhs_m;
//   /// \brief `rhs_m` with Jacobi preconditioner applied.
//   Eigen::Matrix<double, kParameters + kConstraints, 1> rhs_precond_m;
// };

}  // namespace IRL

#include "irl/optimization/constrained_levenberg_marquardt.tpp"

#endif  // IRL_OPTIMIZATION_CONSTRAINED_LEVENBERG_MARQUARDT_H_
