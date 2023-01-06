// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_H_
#define IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_H_

#include <Eigen/Dense>  // Eigen header
#include "irl/parameters/defined_types.h"

namespace IRL {

class ParaboloidGradientLocalZ {
 public:
  static constexpr UnsignedIndex_t NParameters = 1;
  static constexpr bool has_hessian = false;
  ParaboloidGradientLocalZ(void);
  ParaboloidGradientLocalZ(const double a_value);
  ParaboloidGradientLocalZ(const ParaboloidGradientLocalZ& a_gradient);
  ParaboloidGradientLocalZ& operator+=(
      const ParaboloidGradientLocalZ& a_gradient);
  ParaboloidGradientLocalZ& operator=(
      const ParaboloidGradientLocalZ& a_gradient);
  Eigen::Matrix<double, 1, 1>& getGrad(void);
  const Eigen::Matrix<double, 1, 1>& getGrad(void) const;
  void setGrad(const ParaboloidGradientLocalZ& a_rhs);
  double getGradA(void) const;
  double getGradB(void) const;
  double getGradTx(void) const;
  double getGradTy(void) const;
  double getGradTz(void) const;
  double getGradRx(void) const;
  double getGradRy(void) const;
  double getGradRz(void) const;
  void setGradA(const double a_rhs);
  void setGradB(const double a_rhs);
  void setGradTx(const double a_rhs);
  void setGradTy(const double a_rhs);
  void setGradTz(const double a_rhs);
  void setGradRx(const double a_rhs);
  void setGradRy(const double a_rhs);
  void setGradRz(const double a_rhs);
  ~ParaboloidGradientLocalZ(void) = default;

 private:
  // The array onbly contains gradTz
  Eigen::Matrix<double, 1, 1> gradient_m;
};

ParaboloidGradientLocalZ operator*(const ParaboloidGradientLocalZ& a_gradient,
                                   const double a_rhs);
ParaboloidGradientLocalZ operator*(const double a_rhs,
                                   const ParaboloidGradientLocalZ& a_gradient);
ParaboloidGradientLocalZ operator/(const ParaboloidGradientLocalZ& a_gradient,
                                   const double a_rhs);
ParaboloidGradientLocalZ operator+(const ParaboloidGradientLocalZ& a_gradient1,
                                   const ParaboloidGradientLocalZ& a_gradient2);
ParaboloidGradientLocalZ operator-(const ParaboloidGradientLocalZ& a_gradient1,
                                   const ParaboloidGradientLocalZ& a_gradient2);
ParaboloidGradientLocalZ operator-(const ParaboloidGradientLocalZ& a_gradient);

class ParaboloidGradientLocal {
 public:
  static constexpr UnsignedIndex_t NParameters = 8;
  static constexpr bool has_hessian = false;
  ParaboloidGradientLocal(void);
  ParaboloidGradientLocal(const double a_value);
  ParaboloidGradientLocal(const ParaboloidGradientLocal& a_gradient);
  ParaboloidGradientLocal& operator+=(
      const ParaboloidGradientLocal& a_gradient);
  ParaboloidGradientLocal& operator=(const ParaboloidGradientLocal& a_gradient);
  Eigen::Matrix<double, 8, 1>& getGrad(void);
  const Eigen::Matrix<double, 8, 1>& getGrad(void) const;
  void setGrad(const ParaboloidGradientLocal& a_rhs);
  double getGradA(void) const;
  double getGradB(void) const;
  double getGradTx(void) const;
  double getGradTy(void) const;
  double getGradTz(void) const;
  double getGradRx(void) const;
  double getGradRy(void) const;
  double getGradRz(void) const;
  void setGradA(const double a_rhs);
  void setGradB(const double a_rhs);
  void setGradTx(const double a_rhs);
  void setGradTy(const double a_rhs);
  void setGradTz(const double a_rhs);
  void setGradRx(const double a_rhs);
  void setGradRy(const double a_rhs);
  void setGradRz(const double a_rhs);
  ~ParaboloidGradientLocal(void) = default;

 private:
  // The array contains:
  // 0 - gradA
  // 1 - gradB
  // 2 - gradTx
  // 3 - gradTy
  // 4 - gradTz
  // 5 - gradRx
  // 6 - gradRy
  // 7 - gradRz
  Eigen::Matrix<double, 8, 1> gradient_m;
};

ParaboloidGradientLocal operator*(const ParaboloidGradientLocal& a_gradient,
                                  const double a_rhs);
ParaboloidGradientLocal operator*(const double a_rhs,
                                  const ParaboloidGradientLocal& a_gradient);
ParaboloidGradientLocal operator/(const ParaboloidGradientLocal& a_gradient,
                                  const double a_rhs);
ParaboloidGradientLocal operator+(const ParaboloidGradientLocal& a_gradient1,
                                  const ParaboloidGradientLocal& a_gradient2);
ParaboloidGradientLocal operator-(const ParaboloidGradientLocal& a_gradient1,
                                  const ParaboloidGradientLocal& a_gradient2);
ParaboloidGradientLocal operator-(const ParaboloidGradientLocal& a_gradient);

}  // namespace IRL

#include "irl/paraboloid_reconstruction/gradient_paraboloid.tpp"

#endif  // IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_H_
