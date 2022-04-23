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

inline ParaboloidGradientLocalZ::ParaboloidGradientLocalZ(void)
    : gradient_m{0.0} {}

inline ParaboloidGradientLocalZ::ParaboloidGradientLocalZ(const double a_value)
    : gradient_m({a_value}) {}

inline ParaboloidGradientLocalZ ParaboloidGradientLocalZ::operator*(
    const double a_rhs) {
  ParaboloidGradientLocalZ gradient;
  gradient.setGradTz(this->getGradTz() * a_rhs);
  return gradient;
}

inline ParaboloidGradientLocalZ& ParaboloidGradientLocalZ::operator*=(
    const double a_rhs) {
  this->setGradTz(this->getGradTz() * a_rhs);
  return (*this);
}

inline ParaboloidGradientLocalZ ParaboloidGradientLocalZ::operator/(
    const double a_rhs) {
  ParaboloidGradientLocalZ gradient;
  gradient.setGradTz(this->getGradTz() / a_rhs);
  return gradient;
}

inline ParaboloidGradientLocalZ& ParaboloidGradientLocalZ::operator/=(
    const double a_rhs) {
  this->setGradTz(this->getGradTz() / a_rhs);
  return (*this);
}

inline ParaboloidGradientLocalZ ParaboloidGradientLocalZ::operator+(
    const ParaboloidGradientLocalZ& a_gradient) {
  ParaboloidGradientLocalZ gradient;
  const double sumTz = this->getGradTz() + a_gradient.getGradTz();
  gradient.setGradTz(sumTz);
  return gradient;
}

inline ParaboloidGradientLocalZ& ParaboloidGradientLocalZ::operator+=(
    const ParaboloidGradientLocalZ& a_gradient) {
  this->setGradTz(this->getGradTz() + a_gradient.getGradTz());
  return (*this);
}

inline ParaboloidGradientLocalZ ParaboloidGradientLocalZ::operator-(
    const ParaboloidGradientLocalZ& a_gradient) {
  ParaboloidGradientLocalZ gradient;
  gradient.setGradTz(this->getGradTz() - a_gradient.getGradTz());
  return gradient;
}

inline ParaboloidGradientLocalZ& ParaboloidGradientLocalZ::operator-=(
    const ParaboloidGradientLocalZ& a_gradient) {
  this->setGradTz(this->getGradTz() - a_gradient.getGradTz());
  return (*this);
}

inline ParaboloidGradientLocalZ ParaboloidGradientLocalZ::operator=(
    const ParaboloidGradientLocalZ& a_gradient) {
  ParaboloidGradientLocalZ gradient;
  gradient.setGradTz(a_gradient.getGradTz());
  return gradient;
}

inline ParaboloidGradientLocalZ ParaboloidGradientLocalZ::operator=(
    const double a_value) {
  ParaboloidGradientLocalZ gradient;
  gradient.setGradTz(a_value);
  return gradient;
}

inline const double ParaboloidGradientLocalZ::getGradA(void) { return 0.0; }
inline const double ParaboloidGradientLocalZ::getGradB(void) { return 0.0; }
inline const double ParaboloidGradientLocalZ::getGradTx(void) { return 0.0; }
inline const double ParaboloidGradientLocalZ::getGradTy(void) { return 0.0; }
inline const double ParaboloidGradientLocalZ::getGradTz(void) {
  return gradient_m[0];
}
inline const double ParaboloidGradientLocalZ::getGradRx(void) { return 0.0; }
inline const double ParaboloidGradientLocalZ::getGradRy(void) { return 0.0; }
inline const double ParaboloidGradientLocalZ::getGradRz(void) { return 0.0; }
inline const double ParaboloidGradientLocalZ::getGradA(void) const {
  return 0.0;
}
inline const double ParaboloidGradientLocalZ::getGradB(void) const {
  return 0.0;
}
inline const double ParaboloidGradientLocalZ::getGradTx(void) const {
  return 0.0;
}
inline const double ParaboloidGradientLocalZ::getGradTy(void) const {
  return 0.0;
}
inline const double ParaboloidGradientLocalZ::getGradTz(void) const {
  return gradient_m[0];
}
inline const double ParaboloidGradientLocalZ::getGradRx(void) const {
  return 0.0;
}
inline const double ParaboloidGradientLocalZ::getGradRy(void) const {
  return 0.0;
}
inline const double ParaboloidGradientLocalZ::getGradRz(void) const {
  return 0.0;
}
inline void ParaboloidGradientLocalZ::setGradA(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradB(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradTx(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradTy(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradTz(const double a_rhs) {
  gradient_m[0] = a_rhs;
}
inline void ParaboloidGradientLocalZ::setGradRx(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradRy(const double a_rhs) {}
inline void ParaboloidGradientLocalZ::setGradRz(const double a_rhs) {}
inline const bool ParaboloidGradientLocalZ::isLocal(void) const { return true; }

}  // namespace IRL

#endif  // IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_TPP_
