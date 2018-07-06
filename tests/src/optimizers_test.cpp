// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/optimization/optimizers.h"

#include <cmath>

#include "gtest/gtest.h"

#include <Eigen/Dense>  // Eigen header

#include "src/parameters/defined_types.h"

namespace {

using namespace IRL;

// Simple class to look for zeros on with optimization
class PolynomialLM {
  friend class PolynomialOther;

 public:
  PolynomialLM() = default;

  void setup(double a_x) { best_guess_for_x_m(0) = a_x; }

  double calculateY(void) {
    double value = 0.0;
    for (int n = 0; n < order_m + 1; ++n) {
      value += coefficients_m[n] *
               std::pow(guess_for_x_m(0), static_cast<double>(n));
    }
    return value;
  }

  Eigen::Matrix<double, 1, 1> calculateVectorError(void) {
    Eigen::Matrix<double, 1, 1> error_to_return;
    error_to_return(0) = -this->calculateY();
    return error_to_return;
  }

  double calculateScalarError(void) {
    return std::fabs(this->calculateSignedScalarError());
  }

  double calculateSignedScalarError(void) {
    return this->calculateVectorError()(0);
  }

  void updateGuess(Eigen::Matrix<double, 1, 1>* a_update) {
    guess_for_x_m(0) = best_guess_for_x_m(0) + (*a_update)(0);
    current_y_m(0) = this->calculateY();
  }

  void updateBestGuess() {
    best_guess_for_x_m = guess_for_x_m;
    best_y_m = current_y_m;
  }

  bool errorTooHigh(const double a_error) {
    return std::fabs(a_error) > 1.0e-4;
  }

  bool iterationTooHigh(const UnsignedIndex_t a_iteration) {
    return a_iteration > 100;
  }

  bool minimumReached(const Eigen::Matrix<double, 1, 1>& delta_m) {
    return std::fabs(delta_m(0)) < 1.0e-8;
  }

  void increaseLambda(double* a_lambda) { *a_lambda *= lambda_up; }
  void decreaseLambda(double* a_lambda) { *a_lambda *= lambda_down; }

  bool shouldComputeJacobian(const UnsignedIndex_t a_iteration,
                             const UnsignedIndex_t a_last_jacobian) {
    return true;
  }

  Eigen::Matrix<double, 1, 1> calculateChangeInGuess(void) {
    return current_y_m - best_y_m;
  }

  double getSolution(void) { return best_guess_for_x_m(0); }

  ~PolynomialLM() = default;

 private:
  // (x+14)*(x-4)*(x+8)
  int order_m = 3;
  double coefficients_m[4] = {-448.0, 24.0, 18.0, 1.0};
  Eigen::Matrix<double, 1, 1> guess_for_x_m;
  Eigen::Matrix<double, 1, 1> best_guess_for_x_m;
  Eigen::Matrix<double, 1, 1> current_y_m;
  Eigen::Matrix<double, 1, 1> best_y_m;
  static constexpr double lambda_up = 10.0;
  static constexpr double lambda_down = 1.0 / 10.0;
};

class PolynomialOther : public PolynomialLM {
 public:
  PolynomialOther() = default;

  void setup(double a_x) { guess_for_x_m(0) = a_x; }

  void updateGuess(double* a_update) {
    guess_for_x_m(0) = guess_for_x_m(0) + *a_update;
    current_y_m(0) = this->calculateY();
  }

  void setGuess(double* a_guess) {
    guess_for_x_m(0) = *a_guess;
    current_y_m(0) = this->calculateY();
  }

  double getSolution(void) { return guess_for_x_m(0); }

  ~PolynomialOther() = default;
};

TEST(Optimizers, LevenbergMarquardt) {
  PolynomialLM poly;
  poly.setup(0.0);
  Eigen::Matrix<double, 1, 1> initial_delta;
  initial_delta(0) = 0.1;
  LevenbergMarquardt<PolynomialLM, 1, 1> lm_solver;
  lm_solver.solve(&poly, initial_delta);
  EXPECT_NEAR(poly.getSolution(), 4.0, 1.0e-4)
      << "Reason for exit: " << lm_solver.getReason() << '\n';

  poly.setup(-5.0);
  initial_delta(0) = 0.1;
  lm_solver.solve(&poly, initial_delta);
  EXPECT_NEAR(poly.getSolution(), -8.0, 1.0e-4)
      << "Reason for exit: " << lm_solver.getReason() << '\n';

  poly.setup(-20.0);
  initial_delta(0) = 0.1;
  lm_solver.solve(&poly, initial_delta);
  EXPECT_NEAR(poly.getSolution(), -14.0, 1.0e-4)
      << "Reason for exit: " << lm_solver.getReason() << '\n';
}

TEST(Optimizers, BFGS) {
  PolynomialLM poly;
  Eigen::Matrix<double, 1, 1> initial_delta;
  initial_delta(0) = 0.00001;
  BFGS<PolynomialLM, 1> bfgs_solver;

  poly.setup(0.0);
  bfgs_solver.solve(&poly, initial_delta);
  EXPECT_NEAR(poly.getSolution(), 4.0, 1.0e-4)
      << "Reason for exit: " << bfgs_solver.getReason() << '\n';

  poly.setup(-5.0);
  bfgs_solver.solve(&poly, initial_delta);
  EXPECT_NEAR(poly.getSolution(), -8.0, 1.0e-4)
      << "Reason for exit: " << bfgs_solver.getReason() << '\n';

  poly.setup(-50.0);
  bfgs_solver.solve(&poly, initial_delta);
  EXPECT_NEAR(poly.getSolution(), -14.0, 1.0e-4)
      << "Reason for exit: " << bfgs_solver.getReason() << '\n';
}

TEST(Optimizers, Secant) {
  PolynomialOther poly;
  poly.setup(0.0);
  double initial_delta = 0.1;
  Secant<PolynomialOther> sec_solver;
  sec_solver.solve(&poly, initial_delta);
  EXPECT_NEAR(poly.getSolution(), 4.0, 1.0e-4)
      << "Reason for exit: " << sec_solver.getReason() << '\n';

  poly.setup(-5.0);
  initial_delta = 0.1;
  sec_solver.solve(&poly, initial_delta);
  EXPECT_NEAR(poly.getSolution(), -8.0, 1.0e-4)
      << "Reason for exit: " << sec_solver.getReason() << '\n';

  poly.setup(-20.0);
  initial_delta = 0.1;
  sec_solver.solve(&poly, initial_delta);
  EXPECT_NEAR(poly.getSolution(), -14.0, 1.0e-4)
      << "Reason for exit: " << sec_solver.getReason() << '\n';
}

TEST(Optimizers, Illinois) {
  PolynomialOther poly;

  Illinois<PolynomialOther> bi_solver;
  bi_solver.solve(&poly, 0.0, 8.0);
  EXPECT_NEAR(poly.getSolution(), 4.0, 1.0e-4);

  bi_solver.solve(&poly, -13.5, 0.0);
  EXPECT_NEAR(poly.getSolution(), -8.0, 1.0e-4);

  bi_solver.solve(&poly, -30.0, -8.1);
  EXPECT_NEAR(poly.getSolution(), -14.0, 1.0e-4);
}

TEST(Optimizers, Bisection) {
  PolynomialOther poly;

  Bisection<PolynomialOther> bi_solver;
  bi_solver.solve(&poly, 0.0, 8.0);
  EXPECT_NEAR(poly.getSolution(), 4.0, 1.0e-4);

  bi_solver.solve(&poly, -13.5, 0.0);
  EXPECT_NEAR(poly.getSolution(), -8.0, 1.0e-4);

  bi_solver.solve(&poly, -30.0, -8.1);
  EXPECT_NEAR(poly.getSolution(), -14.0, 1.0e-4);
}

TEST(Optimizers, BrentsMethod) {
  PolynomialOther poly;

  BrentsMethod<PolynomialOther> bm_solver;
  bm_solver.solve(&poly, 0.0, 8.0, 1.0e-4);
  EXPECT_NEAR(poly.getSolution(), 4.0, 1.0e-4);

  bm_solver.solve(&poly, -13.5, 0.0, 1.0e-4);
  EXPECT_NEAR(poly.getSolution(), -8.0, 1.0e-4);

  bm_solver.solve(&poly, -30.0, -8.1, 1.0e-4);
  EXPECT_NEAR(poly.getSolution(), -14.0, 1.0e-4);
}

}  // namespace
