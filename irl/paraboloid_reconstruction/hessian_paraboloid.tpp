// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_HESSIAN_PARABOLOID_TPP_
#define IRL_PARABOLOID_RECONSTRUCTION_HESSIAN_PARABOLOID_TPP_

namespace IRL {

template <class GradientType>
inline ParaboloidGradientAndHessianLocal<
    GradientType>::ParaboloidGradientAndHessianLocal(void) {
  gradient_m.setZero();
  hessian_m.setZero();
}

template <class GradientType>
inline ParaboloidGradientAndHessianLocal<
    GradientType>::ParaboloidGradientAndHessianLocal(const double a_value) {
  gradient_m.setConstant(a_value);
  hessian_m.setZero();
}

template <class GradientType>
inline ParaboloidGradientAndHessianLocal<GradientType>::
    ParaboloidGradientAndHessianLocal(
        const ParaboloidGradientAndHessianLocal<GradientType>& a_rhs) {
  gradient_m = a_rhs.getGrad();
  hessian_m = a_rhs.getHessian();
}

template <class GradientType>
inline ParaboloidGradientAndHessianLocal<GradientType>&
ParaboloidGradientAndHessianLocal<GradientType>::operator+=(
    const ParaboloidGradientAndHessianLocal<GradientType>& a_rhs) {
  this->getGrad() += a_rhs.getGrad();
  this->getHessian() += a_rhs.getHessian();
  return (*this);
}

template <class GradientType>
inline ParaboloidGradientAndHessianLocal<GradientType>&
ParaboloidGradientAndHessianLocal<GradientType>::operator=(
    const ParaboloidGradientAndHessianLocal<GradientType>& a_rhs) {
  this->getGrad() = a_rhs.getGrad();
  this->getHessian() = a_rhs.getHessian();
  return (*this);
}

template <class GradientType>
inline ParaboloidGradientAndHessianLocal<GradientType> operator*(
    const ParaboloidGradientAndHessianLocal<GradientType>& a_rhs,
    const double a_factor) {
  ParaboloidGradientAndHessianLocal<GradientType> gradandhess;
  gradandhess.getGrad() = a_rhs.getGrad() * a_factor;
  gradandhess.getHessian() = a_rhs.getHessian() * a_factor;
  return gradandhess;
}

template <class GradientType>
inline ParaboloidGradientAndHessianLocal<GradientType> operator*(
    const double a_factor,
    const ParaboloidGradientAndHessianLocal<GradientType>& a_rhs) {
  return a_rhs * a_factor;
}

template <class GradientType>
inline ParaboloidGradientAndHessianLocal<GradientType> operator/(
    const ParaboloidGradientAndHessianLocal<GradientType>& a_rhs,
    const double a_dividor) {
  ParaboloidGradientAndHessianLocal<GradientType> gradandhess;
  gradandhess.getGrad() = a_rhs.getGrad() / a_dividor;
  gradandhess.getHessian() = a_rhs.getHessian() / a_dividor;
  return gradandhess;
}

template <class GradientType>
inline ParaboloidGradientAndHessianLocal<GradientType> operator+(
    const ParaboloidGradientAndHessianLocal<GradientType>& a_rhs1,
    const ParaboloidGradientAndHessianLocal<GradientType>& a_rhs2) {
  ParaboloidGradientAndHessianLocal<GradientType> gradandhess;
  gradandhess.getGrad() = a_rhs1.getGrad() + a_rhs2.getGrad();
  gradandhess.getHessian() = a_rhs1.getHessian() + a_rhs2.getHessian();
  return gradandhess;
}

template <class GradientType>
inline ParaboloidGradientAndHessianLocal<GradientType> operator-(
    const ParaboloidGradientAndHessianLocal<GradientType>& a_rhs1,
    const ParaboloidGradientAndHessianLocal<GradientType>& a_rhs2) {
  ParaboloidGradientAndHessianLocal<GradientType> gradandhess;
  gradandhess.getGrad() = a_rhs1.getGrad() - a_rhs2.getGrad();
  gradandhess.getHessian() = a_rhs1.getHessian() - a_rhs2.getHessian();
  return gradandhess;
}

template <class GradientType>
inline ParaboloidGradientAndHessianLocal<GradientType> operator-(
    const ParaboloidGradientAndHessianLocal<GradientType>& a_rhs) {
  ParaboloidGradientAndHessianLocal<GradientType> gradandhess;
  gradandhess.getGrad() = -a_rhs.getGrad();
  gradandhess.getHessian() = -a_rhs.getHessian();
  return gradandhess;
}

template <class GradientType>
inline Eigen::Matrix<double, GradientType::NParameters, 1>&
ParaboloidGradientAndHessianLocal<GradientType>::getGrad(void) {
  return gradient_m;
}

template <class GradientType>
inline const Eigen::Matrix<double, GradientType::NParameters, 1>&
ParaboloidGradientAndHessianLocal<GradientType>::getGrad(void) const {
  return gradient_m;
}

template <class GradientType>
inline Eigen::Matrix<double, GradientType::NParameters,
                     GradientType::NParameters>&
ParaboloidGradientAndHessianLocal<GradientType>::getHessian(void) {
  return hessian_m;
}

template <class GradientType>
inline const Eigen::Matrix<double, GradientType::NParameters,
                           GradientType::NParameters>&
ParaboloidGradientAndHessianLocal<GradientType>::getHessian(void) const {
  return hessian_m;
}

template <class GradientType>
inline void ParaboloidGradientAndHessianLocal<GradientType>::setGradA(
    const double a_rhs) {
  gradient_m(0) = a_rhs;
}
template <class GradientType>
inline void ParaboloidGradientAndHessianLocal<GradientType>::setGradB(
    const double a_rhs) {
  gradient_m(1) = a_rhs;
}
template <class GradientType>
inline void ParaboloidGradientAndHessianLocal<GradientType>::setGradTx(
    const double a_rhs) {
  gradient_m(2) = a_rhs;
}
template <class GradientType>
inline void ParaboloidGradientAndHessianLocal<GradientType>::setGradTy(
    const double a_rhs) {
  gradient_m(3) = a_rhs;
}
template <class GradientType>
inline void ParaboloidGradientAndHessianLocal<GradientType>::setGradTz(
    const double a_rhs) {
  gradient_m(4) = a_rhs;
}
template <class GradientType>
inline void ParaboloidGradientAndHessianLocal<GradientType>::setGradRx(
    const double a_rhs) {
  gradient_m(5) = a_rhs;
}
template <class GradientType>
inline void ParaboloidGradientAndHessianLocal<GradientType>::setGradRy(
    const double a_rhs) {
  gradient_m(6) = a_rhs;
}
template <class GradientType>
inline void ParaboloidGradientAndHessianLocal<GradientType>::setGradRz(
    const double a_rhs) {
  gradient_m(7) = a_rhs;
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_TPP_
