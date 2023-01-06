// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GEOMETRY_GENERAL_SCALAR_WITH_GRADIENT_H_
#define IRL_GEOMETRY_GENERAL_SCALAR_WITH_GRADIENT_H_

namespace IRL {

template <class GradientType>
class ScalarWithGradient {
 public:
  using gradient_type = GradientType;
  ScalarWithGradient(void) {
    scalar_m = 0.0;
    gradient_m = GradientType(0.0);
  }
  constexpr ScalarWithGradient(const double a_value) {
    scalar_m = a_value;
    gradient_m = GradientType(0.0);
  }
  constexpr ScalarWithGradient(const double a_value,
                               const GradientType& a_gradient) {
    scalar_m = a_value;
    gradient_m = a_gradient;
  }
  constexpr ScalarWithGradient(
      const ScalarWithGradient<GradientType>& a_scalar) {
    scalar_m = a_scalar.value();
    gradient_m = a_scalar.gradient();
  }
  double& value(void) { return scalar_m; }
  const double& value(void) const { return scalar_m; }
  GradientType& gradient(void) { return gradient_m; }
  const GradientType& gradient(void) const { return gradient_m; }
  // operator double() const { return scalar_m; }
  ScalarWithGradient<GradientType>& operator+=(
      const ScalarWithGradient<GradientType>& a_rhs) {
    scalar_m += a_rhs.scalar_m;
    gradient_m = gradient_m + a_rhs.gradient_m;
    return (*this);
  }
  ScalarWithGradient<GradientType>& operator*=(
      const ScalarWithGradient<GradientType>& a_rhs) {
    scalar_m *= a_rhs.scalar_m;
    gradient_m = gradient_m * a_rhs.scalar_m + scalar_m * a_rhs.gradient_m;
    return (*this);
  }
  ScalarWithGradient<GradientType>& operator*=(const double a_rhs) {
    scalar_m *= a_rhs;
    gradient_m = gradient_m * a_rhs;
    return (*this);
  }
  ~ScalarWithGradient(void) = default;

 private:
  double scalar_m;
  GradientType gradient_m;
};

template <class GradientType>
struct has_embedded_gradient<ScalarWithGradient<GradientType>>
    : std::true_type {};

template <class GradientType>
ScalarWithGradient<GradientType> operator*(
    const ScalarWithGradient<GradientType>& a_scalar, const double a_rhs) {
  ScalarWithGradient<GradientType> new_scalar(a_scalar);
  new_scalar.value() *= a_rhs;
  new_scalar.gradient() = new_scalar.gradient() * a_rhs;
  return new_scalar;
}
template <class GradientType>
ScalarWithGradient<GradientType> operator*(
    const ScalarWithGradient<GradientType>& a_scalar1,
    const ScalarWithGradient<GradientType>& a_scalar2) {
  ScalarWithGradient<GradientType> new_scalar(0.0);
  new_scalar.value() = a_scalar1.value() * a_scalar2.value();
  if constexpr (GradientType::has_hessian) {
    new_scalar.gradient().getGrad() =
        a_scalar1.value() * a_scalar2.gradient().getGrad() +
        a_scalar1.gradient().getGrad() * a_scalar2.value();
    new_scalar.gradient().getHessian() =
        a_scalar1.value() * a_scalar2.gradient().getHessian() +
        a_scalar1.gradient().getHessian() * a_scalar2.value() +
        a_scalar1.gradient().getGrad() *
            a_scalar2.gradient().getGrad().transpose() +
        a_scalar2.gradient().getGrad() *
            a_scalar1.gradient().getGrad().transpose();
  } else {
    new_scalar.gradient() = a_scalar1.value() * a_scalar2.gradient() +
                            a_scalar1.gradient() * a_scalar2.value();
  }
  return new_scalar;
}
template <class GradientType>
ScalarWithGradient<GradientType> operator*(
    const double a_rhs, const ScalarWithGradient<GradientType>& a_scalar) {
  return ScalarWithGradient<GradientType>(a_scalar * a_rhs);
}
template <class GradientType>
ScalarWithGradient<GradientType> operator/(
    const ScalarWithGradient<GradientType>& a_scalar, const double a_rhs) {
  ScalarWithGradient<GradientType> new_scalar(a_scalar);
  new_scalar.value() /= a_rhs;
  new_scalar.gradient() = new_scalar.gradient() / a_rhs;
  return new_scalar;
}
template <class GradientType>
ScalarWithGradient<GradientType> operator/(
    const ScalarWithGradient<GradientType>& a_scalar1,
    const ScalarWithGradient<GradientType>& a_scalar2) {
  ScalarWithGradient<GradientType> new_scalar(0.0);
  new_scalar.value() = a_scalar1.value() / a_scalar2.value();
  if constexpr (GradientType::has_hessian) {
    ScalarWithGradient<GradientType> inv_scalar2(0.0);
    inv_scalar2.value() = 1.0 / a_scalar2.value();
    inv_scalar2.gradient().getGrad() = -a_scalar2.gradient().getGrad() /
                                       (a_scalar2.value() * a_scalar2.value());
    inv_scalar2.gradient().getHessian() =
        -a_scalar2.gradient().getHessian() /
            (a_scalar2.value() * a_scalar2.value()) +
        2.0 * a_scalar2.gradient().getGrad() *
            a_scalar2.gradient().getGrad().transpose() /
            (a_scalar2.value() * a_scalar2.value() * a_scalar2.value());
    new_scalar.gradient().getGrad() =
        a_scalar1.value() * inv_scalar2.gradient().getGrad() +
        a_scalar1.gradient().getGrad() * inv_scalar2.value();
    new_scalar.gradient().getHessian() =
        a_scalar1.value() * inv_scalar2.gradient().getHessian() +
        a_scalar1.gradient().getHessian() * inv_scalar2.value() +
        a_scalar1.gradient().getGrad() *
            inv_scalar2.gradient().getGrad().transpose() +
        inv_scalar2.gradient().getGrad() *
            a_scalar1.gradient().getGrad().transpose();
  } else {
    new_scalar.gradient() = a_scalar1.gradient() / a_scalar2.value() -
                            a_scalar1.value() * a_scalar2.gradient() /
                                (a_scalar2.value() * a_scalar2.value());
  }
  return new_scalar;
}
template <class GradientType>
ScalarWithGradient<GradientType> operator+(
    const ScalarWithGradient<GradientType>& a_scalar1,
    const ScalarWithGradient<GradientType>& a_scalar2) {
  ScalarWithGradient<GradientType> new_scalar(0.0);
  new_scalar.value() = a_scalar1.value() + a_scalar2.value();
  new_scalar.gradient() = a_scalar1.gradient() + a_scalar2.gradient();
  return new_scalar;
}
template <class GradientType>
ScalarWithGradient<GradientType> operator-(
    const ScalarWithGradient<GradientType>& a_scalar1,
    const ScalarWithGradient<GradientType>& a_scalar2) {
  ScalarWithGradient<GradientType> new_scalar(0.0);
  new_scalar.value() = a_scalar1.value() - a_scalar2.value();
  new_scalar.gradient() = a_scalar1.gradient() - a_scalar2.gradient();
  return new_scalar;
}
template <class GradientType>
ScalarWithGradient<GradientType> operator-(
    const ScalarWithGradient<GradientType>& a_scalar) {
  ScalarWithGradient<GradientType> new_scalar(0.0);
  new_scalar.value() = -a_scalar.value();
  new_scalar.gradient() = -a_scalar.gradient();
  return new_scalar;
}
template <class GradientType>
bool operator>=(const ScalarWithGradient<GradientType>& a_scalar1,
                const ScalarWithGradient<GradientType>& a_scalar2) {
  return a_scalar1.value() >= a_scalar2.value();
}
template <class GradientType>
bool operator>(const ScalarWithGradient<GradientType>& a_scalar1,
               const ScalarWithGradient<GradientType>& a_scalar2) {
  return a_scalar1.value() > a_scalar2.value();
}
template <class GradientType>
bool operator<=(const ScalarWithGradient<GradientType>& a_scalar1,
                const ScalarWithGradient<GradientType>& a_scalar2) {
  return a_scalar1.value() <= a_scalar2.value();
}
template <class GradientType>
bool operator<(const ScalarWithGradient<GradientType>& a_scalar1,
               const ScalarWithGradient<GradientType>& a_scalar2) {
  return a_scalar1.value() < a_scalar2.value();
}
template <class ScalarType>
enable_if_t<!has_embedded_gradient<ScalarType>::value, ScalarType> SqrtMoments(
    const ScalarType& a_scalar);
template <class ScalarType>
enable_if_t<!has_embedded_gradient<ScalarType>::value, ScalarType> LogMoments(
    const ScalarType& a_scalar);
template <class ScalarType>
enable_if_t<!has_embedded_gradient<ScalarType>::value, ScalarType>
ArctanMoments(const ScalarType& a_scalar);
template <class ScalarType>
enable_if_t<!has_embedded_gradient<ScalarType>::value, ScalarType>
ArctanhMoments(const ScalarType& a_scalar);
template <class ScalarType>
enable_if_t<!has_embedded_gradient<ScalarType>::value, ScalarType> PowMoments(
    const ScalarType& a_scalar, const double a_power);
template <class ScalarType>
inline enable_if_t<has_embedded_gradient<ScalarType>::value, ScalarType>
SqrtMoments(const ScalarType& a_scalar) {
  ScalarType new_scalar(0.0);
  new_scalar.value() = std::sqrt(a_scalar.value());
  if constexpr (ScalarType::gradient_type::has_hessian) {
    new_scalar.gradient().getGrad() =
        a_scalar.gradient().getGrad() / (2.0 * new_scalar.value());
    new_scalar.gradient().getHessian() =
        a_scalar.gradient().getHessian() / (2.0 * new_scalar.value()) -
        a_scalar.gradient().getGrad() *
            new_scalar.gradient().getGrad().transpose() /
            (2.0 * new_scalar.value() * new_scalar.value());
  } else {
    new_scalar.gradient() = a_scalar.gradient() / (2.0 * new_scalar.value());
  }
  return new_scalar;
}
template <>
inline enable_if_t<!has_embedded_gradient<double>::value, double> SqrtMoments(
    const double& a_scalar) {
  return std::sqrt(a_scalar);
}
template <class ScalarType>
inline enable_if_t<has_embedded_gradient<ScalarType>::value, ScalarType>
LogMoments(const ScalarType& a_scalar) {
  ScalarType new_scalar(0.0);
  new_scalar.value() = std::log(a_scalar.value());
  if constexpr (ScalarType::gradient_type::has_hessian) {
    new_scalar.gradient().getGrad() =
        a_scalar.gradient().getGrad() / a_scalar.value();
    new_scalar.gradient().getHessian() =
        a_scalar.gradient().getHessian() / a_scalar.value() -
        a_scalar.gradient().getGrad() *
            a_scalar.gradient().getGrad().transpose() /
            (a_scalar.value() * a_scalar.value());
  } else {
    new_scalar.gradient() = a_scalar.gradient() / a_scalar.value();
  }
  return new_scalar;
}
template <>
inline enable_if_t<!has_embedded_gradient<double>::value, double> LogMoments(
    const double& a_scalar) {
  return std::log(a_scalar);
}
template <class ScalarType>
inline enable_if_t<has_embedded_gradient<ScalarType>::value, ScalarType>
ArctanMoments(const ScalarType& a_scalar) {
  ScalarType new_scalar(0.0);
  new_scalar.value() = std::atan(a_scalar.value());
  if constexpr (ScalarType::gradient_type::has_hessian) {
    new_scalar.gradient().getGrad() =
        a_scalar.gradient().getGrad() /
        (1.0 + a_scalar.value() * a_scalar.value());
    new_scalar.gradient().getHessian() =
        a_scalar.gradient().getHessian() /
            (1.0 + a_scalar.value() * a_scalar.value()) -
        a_scalar.gradient().getGrad() *
            (2.0 * a_scalar.gradient().getGrad().transpose() *
             a_scalar.value()) /
            ((1.0 + a_scalar.value() * a_scalar.value()) *
             (1.0 + a_scalar.value() * a_scalar.value()));
  } else {
    new_scalar.gradient() =
        a_scalar.gradient() / (1.0 + a_scalar.value() * a_scalar.value());
  }
  return new_scalar;
}
template <>
inline enable_if_t<!has_embedded_gradient<double>::value, double> ArctanMoments(
    const double& a_scalar) {
  return std::atan(a_scalar);
}
template <class ScalarType>
inline enable_if_t<has_embedded_gradient<ScalarType>::value, ScalarType>
ArctanhMoments(const ScalarType& a_scalar) {
  ScalarType new_scalar(0.0);
  new_scalar.value() = std::atanh(a_scalar.value());
  if constexpr (ScalarType::gradient_type::has_hessian) {
    new_scalar.gradient().getGrad() =
        a_scalar.gradient().getGrad() /
        (1.0 - a_scalar.value() * a_scalar.value());
    new_scalar.gradient().getHessian() =
        a_scalar.gradient().getHessian() /
            (1.0 - a_scalar.value() * a_scalar.value()) +
        a_scalar.gradient().getGrad() *
            (2.0 * a_scalar.gradient().getGrad().transpose() *
             a_scalar.value()) /
            ((1.0 - a_scalar.value() * a_scalar.value()) *
             (1.0 - a_scalar.value() * a_scalar.value()));
  } else {
    new_scalar.gradient() =
        a_scalar.gradient() / (1.0 - a_scalar.value() * a_scalar.value());
  }
  return new_scalar;
}
template <>
inline enable_if_t<!has_embedded_gradient<double>::value, double>
ArctanhMoments(const double& a_scalar) {
  return std::atanh(a_scalar);
}
template <class ScalarType>
inline enable_if_t<has_embedded_gradient<ScalarType>::value, ScalarType>
PowMoments(const ScalarType& a_scalar, const double a_power) {
  ScalarType new_scalar(0.0);
  new_scalar.value() = std::pow(a_scalar.value(), a_power);
  if constexpr (ScalarType::gradient_type::has_hessian) {
    new_scalar.gradient().getGrad() = a_power * a_scalar.gradient().getGrad() *
                                      std::pow(a_scalar.value(), a_power - 1.0);
    new_scalar.gradient().getHessian() =
        a_power * a_scalar.gradient().getHessian() *
            std::pow(a_scalar.value(), a_power - 1.0) +
        a_power * a_scalar.gradient().getGrad() * (a_power - 1.0) *
            a_scalar.gradient().getGrad().transpose() *
            std::pow(a_scalar.value(), a_power - 2.0);
  } else {
    new_scalar.gradient() = a_power * a_scalar.gradient() *
                            std::pow(a_scalar.value(), a_power - 1.0);
  }
  return new_scalar;
}
template <>
inline enable_if_t<!has_embedded_gradient<double>::value, double> PowMoments(
    const double& a_scalar, const double a_power) {
  return std::pow(a_scalar, a_power);
}

}  // namespace IRL

#include "irl/geometry/general/scalar_with_gradient.tpp"

#endif