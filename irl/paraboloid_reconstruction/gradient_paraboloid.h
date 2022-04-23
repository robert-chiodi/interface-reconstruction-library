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

#include "irl/parameters/defined_types.h"

namespace IRL {

class ParaboloidGradientLocalZ {
 public:
  ParaboloidGradientLocalZ(void);
  ParaboloidGradientLocalZ(const double a_value);
  inline double& operator[](const UnsignedIndex_t a_elem);
  inline ParaboloidGradientLocalZ operator*(const double a_rhs);
  inline ParaboloidGradientLocalZ& operator*=(const double a_rhs);
  inline ParaboloidGradientLocalZ operator/(const double a_rhs);
  inline ParaboloidGradientLocalZ& operator/=(const double a_rhs);
  inline ParaboloidGradientLocalZ operator+(
      const ParaboloidGradientLocalZ& a_gradient);
  inline ParaboloidGradientLocalZ& operator+=(
      const ParaboloidGradientLocalZ& a_gradient);
  inline ParaboloidGradientLocalZ operator-(
      const ParaboloidGradientLocalZ& a_gradient);
  inline ParaboloidGradientLocalZ& operator-=(
      const ParaboloidGradientLocalZ& a_gradient);
  inline ParaboloidGradientLocalZ operator=(
      const ParaboloidGradientLocalZ& a_gradient);
  inline ParaboloidGradientLocalZ operator=(const double a_value);
  inline const double getGradA(void);
  inline const double getGradB(void);
  inline const double getGradTx(void);
  inline const double getGradTy(void);
  inline const double getGradTz(void);
  inline const double getGradRx(void);
  inline const double getGradRy(void);
  inline const double getGradRz(void);
  inline const double getGradA(void) const;
  inline const double getGradB(void) const;
  inline const double getGradTx(void) const;
  inline const double getGradTy(void) const;
  inline const double getGradTz(void) const;
  inline const double getGradRx(void) const;
  inline const double getGradRy(void) const;
  inline const double getGradRz(void) const;
  inline void setGradA(const double a_rhs);
  inline void setGradB(const double a_rhs);
  inline void setGradTx(const double a_rhs);
  inline void setGradTy(const double a_rhs);
  inline void setGradTz(const double a_rhs);
  inline void setGradRx(const double a_rhs);
  inline void setGradRy(const double a_rhs);
  inline void setGradRz(const double a_rhs);
  inline const bool isLocal(void) const;
  ~ParaboloidGradientLocalZ(void) = default;

 private:
  std::array<double, 1> gradient_m;
};

}  // namespace IRL

#include "irl/paraboloid_reconstruction/gradient_paraboloid.tpp"

#endif  // IRL_PARABOLOID_RECONSTRUCTION_GRADIENT_PARABOLOID_H_
