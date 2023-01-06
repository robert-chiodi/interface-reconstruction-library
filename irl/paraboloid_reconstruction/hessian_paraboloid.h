// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_HESSIAN_PARABOLOID_H_
#define IRL_PARABOLOID_RECONSTRUCTION_HESSIAN_PARABOLOID_H_

#include <Eigen/Dense>  // Eigen header
#include "irl/parameters/defined_types.h"

namespace IRL {

template <class GradientType>
class ParaboloidGradientAndHessianLocal {
 public:
  using gradient_type = GradientType;
  static constexpr UnsignedIndex_t NParameters = GradientType::NParameters;
  static constexpr bool has_hessian = true;
  ParaboloidGradientAndHessianLocal(void);
  ParaboloidGradientAndHessianLocal(const double a_value);
  ParaboloidGradientAndHessianLocal(
      const ParaboloidGradientAndHessianLocal<GradientType>& a_gradient);
  ParaboloidGradientAndHessianLocal& operator+=(
      const ParaboloidGradientAndHessianLocal<GradientType>& a_gradient);
  ParaboloidGradientAndHessianLocal& operator=(
      const ParaboloidGradientAndHessianLocal<GradientType>& a_gradient);
  Eigen::Matrix<double, GradientType::NParameters, 1>& getGrad(void);
  const Eigen::Matrix<double, GradientType::NParameters, 1>& getGrad(
      void) const;
  Eigen::Matrix<double, GradientType::NParameters, GradientType::NParameters>&
  getHessian(void);
  const Eigen::Matrix<double, GradientType::NParameters,
                      GradientType::NParameters>&
  getHessian(void) const;
  void setGradA(const double a_rhs);
  void setGradB(const double a_rhs);
  void setGradTx(const double a_rhs);
  void setGradTy(const double a_rhs);
  void setGradTz(const double a_rhs);
  void setGradRx(const double a_rhs);
  void setGradRy(const double a_rhs);
  void setGradRz(const double a_rhs);
  ~ParaboloidGradientAndHessianLocal(void) = default;

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
  Eigen::Matrix<double, NParameters, 1> gradient_m;
  Eigen::Matrix<double, NParameters, NParameters> hessian_m;
};

template <class GradientType>
ParaboloidGradientAndHessianLocal<GradientType> operator*(
    const ParaboloidGradientAndHessianLocal<GradientType>& a_gradient,
    const double a_rhs);
template <class GradientType>
ParaboloidGradientAndHessianLocal<GradientType> operator*(
    const double a_rhs,
    const ParaboloidGradientAndHessianLocal<GradientType>& a_gradient);
template <class GradientType>
ParaboloidGradientAndHessianLocal<GradientType> operator/(
    const ParaboloidGradientAndHessianLocal<GradientType>& a_gradient,
    const double a_rhs);
template <class GradientType>
ParaboloidGradientAndHessianLocal<GradientType> operator+(
    const ParaboloidGradientAndHessianLocal<GradientType>& a_gradient1,
    const ParaboloidGradientAndHessianLocal<GradientType>& a_gradient2);
template <class GradientType>
ParaboloidGradientAndHessianLocal<GradientType> operator-(
    const ParaboloidGradientAndHessianLocal<GradientType>& a_gradient1,
    const ParaboloidGradientAndHessianLocal<GradientType>& a_gradient2);
template <class GradientType>
ParaboloidGradientAndHessianLocal<GradientType> operator-(
    const ParaboloidGradientAndHessianLocal<GradientType>& a_gradient);

}  // namespace IRL

#include "irl/paraboloid_reconstruction/hessian_paraboloid.tpp"

#endif  // IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_H_
