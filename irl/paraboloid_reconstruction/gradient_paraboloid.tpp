// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_TPP_
#define IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_TPP_

namespace IRL {

inline ParaboloidGradientLocalZ::ParaboloidGradientLocalZ(void) {
  gradient_m.setZero();
}

inline ParaboloidGradientLocalZ::ParaboloidGradientLocalZ(
    const double a_value) {
  gradient_m.setConstant(a_value);
}

inline ParaboloidGradientLocalZ::ParaboloidGradientLocalZ(
    const ParaboloidGradientLocalZ& a_gradient) {
  gradient_m = a_gradient.getGrad();
}

inline ParaboloidGradientLocalZ& ParaboloidGradientLocalZ::operator+=(
    const ParaboloidGradientLocalZ& a_gradient) {
  this->getGrad() += a_gradient.getGrad();
  return (*this);
}

inline ParaboloidGradientLocalZ& ParaboloidGradientLocalZ::operator=(
    const ParaboloidGradientLocalZ& a_gradient) {
  this->getGrad() = a_gradient.getGrad();
  return (*this);
}

inline ParaboloidGradientLocalZ operator*(
    const ParaboloidGradientLocalZ& a_gradient, const double a_rhs) {
  ParaboloidGradientLocalZ gradient(0.0);
  gradient.getGrad() = a_gradient.getGrad() * a_rhs;
  return gradient;
}

inline ParaboloidGradientLocalZ operator*(
    const double a_rhs, const ParaboloidGradientLocalZ& a_gradient) {
  return a_gradient * a_rhs;
}

inline ParaboloidGradientLocalZ operator/(
    const ParaboloidGradientLocalZ& a_gradient, const double a_rhs) {
  ParaboloidGradientLocalZ gradient(0.0);
  gradient.getGrad() = a_gradient.getGrad() / a_rhs;
  return gradient;
}

inline ParaboloidGradientLocalZ operator+(
    const ParaboloidGradientLocalZ& a_gradient1,
    const ParaboloidGradientLocalZ& a_gradient2) {
  ParaboloidGradientLocalZ gradient(0.0);
  gradient.getGrad() = a_gradient1.getGrad() + a_gradient2.getGrad();
  return gradient;
}

inline ParaboloidGradientLocalZ operator-(
    const ParaboloidGradientLocalZ& a_gradient1,
    const ParaboloidGradientLocalZ& a_gradient2) {
  ParaboloidGradientLocalZ gradient(0.0);
  gradient.getGrad() = a_gradient1.getGrad() - a_gradient2.getGrad();
  return gradient;
}
inline ParaboloidGradientLocalZ operator-(
    const ParaboloidGradientLocalZ& a_gradient) {
  ParaboloidGradientLocalZ gradient(0.0);
  gradient.getGrad() = -a_gradient.getGrad();
  return gradient;
}

inline Eigen::Matrix<double, 1, 1>& ParaboloidGradientLocalZ::getGrad(void) {
  return gradient_m;
}

inline const Eigen::Matrix<double, 1, 1>& ParaboloidGradientLocalZ::getGrad(
    void) const {
  return gradient_m;
}

inline double ParaboloidGradientLocalZ::getGradA(void) const { return 0.0; }
inline double ParaboloidGradientLocalZ::getGradB(void) const { return 0.0; }
inline double ParaboloidGradientLocalZ::getGradTx(void) const { return 0.0; }
inline double ParaboloidGradientLocalZ::getGradTy(void) const { return 0.0; }
inline double ParaboloidGradientLocalZ::getGradTz(void) const {
  return gradient_m(0);
}
inline double ParaboloidGradientLocalZ::getGradRx(void) const { return 0.0; }
inline double ParaboloidGradientLocalZ::getGradRy(void) const { return 0.0; }
inline double ParaboloidGradientLocalZ::getGradRz(void) const { return 0.0; }
inline void ParaboloidGradientLocalZ::setGrad(
    const ParaboloidGradientLocalZ& a_rhs) {
  gradient_m = a_rhs.getGrad();
}
inline void ParaboloidGradientLocalZ::setGradA(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradB(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradTx(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradTy(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradTz(const double a_rhs) {
  gradient_m(0) = a_rhs;
}
inline void ParaboloidGradientLocalZ::setGradRx(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradRy(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradRz(const double a_rhs) {}

inline ParaboloidGradientLocal::ParaboloidGradientLocal(void) {
  gradient_m.setZero();
}

inline ParaboloidGradientLocal::ParaboloidGradientLocal(const double a_value) {
  gradient_m.setConstant(a_value);
}
inline ParaboloidGradientLocal::ParaboloidGradientLocal(
    const ParaboloidGradientLocal& a_gradient) {
  gradient_m = a_gradient.getGrad();
}

inline ParaboloidGradientLocal& ParaboloidGradientLocal::operator+=(
    const ParaboloidGradientLocal& a_gradient) {
  this->getGrad() += a_gradient.getGrad();
  return (*this);
}

inline ParaboloidGradientLocal& ParaboloidGradientLocal::operator=(
    const ParaboloidGradientLocal& a_gradient) {
  this->getGrad() = a_gradient.getGrad();
  return (*this);
}

inline ParaboloidGradientLocal operator*(
    const ParaboloidGradientLocal& a_gradient, const double a_rhs) {
  ParaboloidGradientLocal gradient(0.0);
  gradient.getGrad() = a_gradient.getGrad() * a_rhs;
  return gradient;
}

inline ParaboloidGradientLocal operator*(
    const double a_rhs, const ParaboloidGradientLocal& a_gradient) {
  return a_gradient * a_rhs;
}

inline ParaboloidGradientLocal operator/(
    const ParaboloidGradientLocal& a_gradient, const double a_rhs) {
  ParaboloidGradientLocal gradient(0.0);
  gradient.getGrad() = a_gradient.getGrad() / a_rhs;
  return gradient;
}

inline ParaboloidGradientLocal operator+(
    const ParaboloidGradientLocal& a_gradient1,
    const ParaboloidGradientLocal& a_gradient2) {
  ParaboloidGradientLocal gradient(0.0);
  gradient.getGrad() = a_gradient1.getGrad() + a_gradient2.getGrad();
  return gradient;
}

inline ParaboloidGradientLocal operator-(
    const ParaboloidGradientLocal& a_gradient1,
    const ParaboloidGradientLocal& a_gradient2) {
  ParaboloidGradientLocal gradient(0.0);
  gradient.getGrad() = a_gradient1.getGrad() - a_gradient2.getGrad();
  return gradient;
}
inline ParaboloidGradientLocal operator-(
    const ParaboloidGradientLocal& a_gradient) {
  ParaboloidGradientLocal gradient(0.0);
  gradient.getGrad() = -a_gradient.getGrad();
  return gradient;
}

inline Eigen::Matrix<double, 8, 1>& ParaboloidGradientLocal::getGrad(void) {
  return gradient_m;
}

inline const Eigen::Matrix<double, 8, 1>& ParaboloidGradientLocal::getGrad(
    void) const {
  return gradient_m;
}

inline double ParaboloidGradientLocal::getGradA(void) const {
  return gradient_m(0);
}
inline double ParaboloidGradientLocal::getGradB(void) const {
  return gradient_m(1);
}
inline double ParaboloidGradientLocal::getGradTx(void) const {
  return gradient_m(2);
}
inline double ParaboloidGradientLocal::getGradTy(void) const {
  return gradient_m(3);
}
inline double ParaboloidGradientLocal::getGradTz(void) const {
  return gradient_m(4);
}
inline double ParaboloidGradientLocal::getGradRx(void) const {
  return gradient_m(5);
}
inline double ParaboloidGradientLocal::getGradRy(void) const {
  return gradient_m(6);
}
inline double ParaboloidGradientLocal::getGradRz(void) const {
  return gradient_m(7);
}
inline void ParaboloidGradientLocal::setGrad(
    const ParaboloidGradientLocal& a_rhs) {
  gradient_m = a_rhs.getGrad();
}
inline void ParaboloidGradientLocal::setGradA(const double a_rhs) {
  gradient_m(0) = a_rhs;
}
inline void ParaboloidGradientLocal::setGradB(const double a_rhs) {
  gradient_m(1) = a_rhs;
}
inline void ParaboloidGradientLocal::setGradTx(const double a_rhs) {
  gradient_m(2) = a_rhs;
}
inline void ParaboloidGradientLocal::setGradTy(const double a_rhs) {
  gradient_m(3) = a_rhs;
}
inline void ParaboloidGradientLocal::setGradTz(const double a_rhs) {
  gradient_m(4) = a_rhs;
}
inline void ParaboloidGradientLocal::setGradRx(const double a_rhs) {
  gradient_m(5) = a_rhs;
}
inline void ParaboloidGradientLocal::setGradRy(const double a_rhs) {
  gradient_m(6) = a_rhs;
}
inline void ParaboloidGradientLocal::setGradRz(const double a_rhs) {
  gradient_m(7) = a_rhs;
}

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_TPP_
