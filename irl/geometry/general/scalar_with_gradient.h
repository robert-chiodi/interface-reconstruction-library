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

template <class ScalarType, class GradientType>
class ScalarWithGradientBase {
 public:
  using value_type = ScalarType;
  using gradient_type = GradientType;
  ScalarWithGradientBase(void) {
    scalar_m = static_cast<ScalarType>(0);
    gradient_m = GradientType(static_cast<ScalarType>(0));
  }
  constexpr ScalarWithGradientBase(const ScalarType a_value) {
    scalar_m = a_value;
    gradient_m = GradientType(static_cast<ScalarType>(0));
  }
  constexpr ScalarWithGradientBase(const ScalarType a_value,
                                   const GradientType& a_gradient) {
    scalar_m = a_value;
    gradient_m = a_gradient;
  }
  constexpr ScalarWithGradientBase(
      const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar) {
    scalar_m = a_scalar.value();
    gradient_m = a_scalar.gradient();
  }
  ScalarType& value(void) { return scalar_m; }
  const ScalarType& value(void) const { return scalar_m; }
  GradientType& gradient(void) { return gradient_m; }
  const GradientType& gradient(void) const { return gradient_m; }
  // operator double() const { return scalar_m; }
  ScalarWithGradientBase<ScalarType, GradientType>& operator+=(
      const ScalarWithGradientBase<ScalarType, GradientType>& a_rhs) {
    scalar_m += a_rhs.scalar_m;
    gradient_m = gradient_m + a_rhs.gradient_m;
    return (*this);
  }
  ScalarWithGradientBase<ScalarType, GradientType>& operator*=(
      const ScalarWithGradientBase<ScalarType, GradientType>& a_rhs) {
    scalar_m *= a_rhs.scalar_m;
    gradient_m = gradient_m * a_rhs.scalar_m + scalar_m * a_rhs.gradient_m;
    return (*this);
  }
  ScalarWithGradientBase<ScalarType, GradientType>& operator*=(
      const ScalarType a_rhs) {
    scalar_m *= a_rhs;
    gradient_m = gradient_m * a_rhs;
    return (*this);
  }
  ~ScalarWithGradientBase(void) = default;

 private:
  ScalarType scalar_m;
  GradientType gradient_m;
};

// template <class ScalarType, class GradientType>
// struct has_embedded_gradient<ScalarWithGradientBase<ScalarType,
// GradientType>>
//     : std::true_type {};

template <class GradientType>
struct has_embedded_gradient<ScalarWithGradientBase<double, GradientType>>
    : std::true_type {};

template <class GradientType>
struct has_embedded_gradient<ScalarWithGradientBase<Quad_t, GradientType>>
    : std::true_type {};

template <class ScalarType, class GradientType>
ScalarWithGradientBase<ScalarType, GradientType> operator*(
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar,
    const ScalarType a_rhs) {
  ScalarWithGradientBase<ScalarType, GradientType> new_scalar(a_scalar);
  new_scalar.value() *= a_rhs;
  new_scalar.gradient() = new_scalar.gradient() * a_rhs;
  return new_scalar;
}

template <class ScalarType, class GradientType>
ScalarWithGradientBase<ScalarType, GradientType> operator*(
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar2) {
  ScalarWithGradientBase<ScalarType, GradientType> new_scalar(
      static_cast<ScalarType>(0));
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

template <class ScalarType, class GradientType>
ScalarWithGradientBase<ScalarType, GradientType> operator*(
    const ScalarType a_rhs,
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar) {
  return ScalarWithGradientBase<ScalarType, GradientType>(a_scalar * a_rhs);
}

template <class ScalarType, class GradientType>
ScalarWithGradientBase<ScalarType, GradientType> operator/(
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar,
    const ScalarType a_rhs) {
  ScalarWithGradientBase<ScalarType, GradientType> new_scalar(a_scalar);
  new_scalar.value() /= a_rhs;
  new_scalar.gradient() = new_scalar.gradient() / a_rhs;
  return new_scalar;
}

template <class ScalarType, class GradientType>
ScalarWithGradientBase<ScalarType, GradientType> operator/(
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar2) {
  ScalarWithGradientBase<ScalarType, GradientType> new_scalar(
      static_cast<ScalarType>(0));
  new_scalar.value() = a_scalar1.value() / a_scalar2.value();
  if constexpr (GradientType::has_hessian) {
    ScalarWithGradientBase<ScalarType, GradientType> inv_scalar2(
        static_cast<ScalarType>(0));
    inv_scalar2.value() = static_cast<ScalarType>(1) / a_scalar2.value();
    inv_scalar2.gradient().getGrad() = -a_scalar2.gradient().getGrad() /
                                       (a_scalar2.value() * a_scalar2.value());
    inv_scalar2.gradient().getHessian() =
        -a_scalar2.gradient().getHessian() /
            (a_scalar2.value() * a_scalar2.value()) +
        static_cast<ScalarType>(2) * a_scalar2.gradient().getGrad() *
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

template <class ScalarType, class GradientType>
ScalarWithGradientBase<ScalarType, GradientType> operator+(
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar2) {
  ScalarWithGradientBase<ScalarType, GradientType> new_scalar(
      static_cast<ScalarType>(0));
  new_scalar.value() = a_scalar1.value() + a_scalar2.value();
  new_scalar.gradient() = a_scalar1.gradient() + a_scalar2.gradient();
  return new_scalar;
}

template <class ScalarType, class GradientType>
ScalarWithGradientBase<ScalarType, GradientType> operator-(
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar2) {
  ScalarWithGradientBase<ScalarType, GradientType> new_scalar(
      static_cast<ScalarType>(0));
  new_scalar.value() = a_scalar1.value() - a_scalar2.value();
  new_scalar.gradient() = a_scalar1.gradient() - a_scalar2.gradient();
  return new_scalar;
}

template <class ScalarType, class GradientType>
ScalarWithGradientBase<ScalarType, GradientType> operator-(
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar) {
  ScalarWithGradientBase<ScalarType, GradientType> new_scalar(
      static_cast<ScalarType>(0));
  new_scalar.value() = -a_scalar.value();
  new_scalar.gradient() = -a_scalar.gradient();
  return new_scalar;
}

template <class ScalarType, class GradientType>
bool operator>=(
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar2) {
  return a_scalar1.value() >= a_scalar2.value();
}

template <class ScalarType, class GradientType>
bool operator>(
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar2) {
  return a_scalar1.value() > a_scalar2.value();
}

template <class ScalarType, class GradientType>
bool operator<=(
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar2) {
  return a_scalar1.value() <= a_scalar2.value();
}

template <class ScalarType, class GradientType>
bool operator<(
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar1,
    const ScalarWithGradientBase<ScalarType, GradientType>& a_scalar2) {
  return a_scalar1.value() < a_scalar2.value();
}

/************************************ SQRT ************************************/

template <class ScalarWithGradientType>
enable_if_t<!has_embedded_gradient<ScalarWithGradientType>::value,
            ScalarWithGradientType>
SqrtMoments(const ScalarWithGradientType& a_scalar) {
  return sqrt(a_scalar);
}

// template <>
// inline enable_if_t<!has_embedded_gradient<double>::value, double>
// SqrtMoments(
//     const double& a_scalar) {
//   return std::sqrt(a_scalar);
// }

// template <>
// inline enable_if_t<!has_embedded_gradient<Quad_t>::value, Quad_t>
// SqrtMoments(
//     const double& a_scalar) {
//   return sqrtq(a_scalar);
// }

template <class ScalarWithGradientType>
inline enable_if_t<has_embedded_gradient<ScalarWithGradientType>::value,
                   ScalarWithGradientType>
SqrtMoments(const ScalarWithGradientType& a_scalar) {
  using ScalarType = typename ScalarWithGradientType::value_type;
  ScalarWithGradientType new_scalar(static_cast<ScalarType>(0));
  new_scalar.value() = sqrt(a_scalar.value());
  if constexpr (ScalarWithGradientType::gradient_type::has_hessian) {
    new_scalar.gradient().getGrad() =
        a_scalar.gradient().getGrad() /
        (static_cast<ScalarType>(2) * new_scalar.value());
    new_scalar.gradient().getHessian() =
        a_scalar.gradient().getHessian() /
            (static_cast<ScalarType>(2) * new_scalar.value()) -
        a_scalar.gradient().getGrad() *
            new_scalar.gradient().getGrad().transpose() /
            (static_cast<ScalarType>(2) * new_scalar.value() *
             new_scalar.value());
  } else {
    new_scalar.gradient() =
        a_scalar.gradient() / (static_cast<ScalarType>(2) * new_scalar.value());
  }
  return new_scalar;
}

/************************************ POW ************************************/

template <class ScalarWithGradientType>
enable_if_t<!has_embedded_gradient<ScalarWithGradientType>::value,
            ScalarWithGradientType>
PowMoments(const ScalarWithGradientType& a_scalar,
           const typename ScalarWithGradientType::value_type a_power) {
  return pow(a_scalar, a_power);
}

// template <>
// inline enable_if_t<!has_embedded_gradient<double>::value, double> PowMoments(
//     const double& a_scalar, const double a_power) {
//   return std::pow(a_scalar, a_power);
// }
// template <>
// inline enable_if_t<!has_embedded_gradient<Quad_t>::value, Quad_t> PowMoments(
//     const Quad_t& a_scalar, const Quad_t a_power) {
//   return powq(a_scalar, a_power);
// }

template <class ScalarWithGradientType>
inline enable_if_t<has_embedded_gradient<ScalarWithGradientType>::value,
                   ScalarWithGradientType>
PowMoments(const ScalarWithGradientType& a_scalar,
           const typename ScalarWithGradientType::value_type a_power) {
  using ScalarType = typename ScalarWithGradientType::value_type;
  ScalarWithGradientType new_scalar(static_cast<ScalarType>(0));
  new_scalar.value() = pow(a_scalar.value(), a_power);
  if constexpr (ScalarWithGradientType::gradient_type::has_hessian) {
    new_scalar.gradient().getGrad() =
        a_power * a_scalar.gradient().getGrad() *
        pow(a_scalar.value(), a_power - static_cast<ScalarType>(1));
    new_scalar.gradient().getHessian() =
        a_power * a_scalar.gradient().getHessian() *
            pow(a_scalar.value(), a_power - static_cast<ScalarType>(1)) +
        a_power * a_scalar.gradient().getGrad() *
            (a_power - static_cast<ScalarType>(1)) *
            a_scalar.gradient().getGrad().transpose() *
            pow(a_scalar.value(), a_power - static_cast<ScalarType>(2));
  } else {
    new_scalar.gradient() =
        a_power * a_scalar.gradient() *
        pow(a_scalar.value(), a_power - static_cast<ScalarType>(1));
  }
  return new_scalar;
}

/************************************ LOG ************************************/

template <class ScalarWithGradientType>
enable_if_t<!has_embedded_gradient<ScalarWithGradientType>::value,
            ScalarWithGradientType>
LogMoments(const ScalarWithGradientType& a_scalar) {
  return log(a_scalar);
}

// template <>
// inline enable_if_t<!has_embedded_gradient<double>::value, double> LogMoments(
//     const double& a_scalar) {
//   return std::log(a_scalar);
// }
// template <>
// inline enable_if_t<!has_embedded_gradient<Quad_t>::value, Quad_t> LogMoments(
//     const Quad_t& a_scalar) {
//   return logq(a_scalar);
// }

template <class ScalarWithGradientType>
inline enable_if_t<has_embedded_gradient<ScalarWithGradientType>::value,
                   ScalarWithGradientType>
LogMoments(const ScalarWithGradientType& a_scalar) {
  using ScalarType = typename ScalarWithGradientType::value_type;
  ScalarWithGradientType new_scalar(static_cast<ScalarType>(0));
  new_scalar.value() = log(a_scalar.value());
  if constexpr (ScalarWithGradientType::gradient_type::has_hessian) {
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

/************************************ ATAN ************************************/

template <class ScalarWithGradientType>
enable_if_t<!has_embedded_gradient<ScalarWithGradientType>::value,
            ScalarWithGradientType>
ArctanMoments(const ScalarWithGradientType& a_scalar) {
  return atan(a_scalar);
}

// template <>
// inline enable_if_t<!has_embedded_gradient<double>::value, double>
// ArctanMoments(
//     const double& a_scalar) {
//   return std::atan(a_scalar);
// }
// template <>
// inline enable_if_t<!has_embedded_gradient<Quad_t>::value, Quad_t>
// ArctanMoments(
//     const Quad_t& a_scalar) {
//   return atanq(a_scalar);
// }

template <class ScalarWithGradientType>
inline enable_if_t<has_embedded_gradient<ScalarWithGradientType>::value,
                   ScalarWithGradientType>
ArctanMoments(const ScalarWithGradientType& a_scalar) {
  using ScalarType = typename ScalarWithGradientType::value_type;
  ScalarWithGradientType new_scalar(static_cast<ScalarType>(0));
  new_scalar.value() = atan(a_scalar.value());
  if constexpr (ScalarWithGradientType::gradient_type::has_hessian) {
    new_scalar.gradient().getGrad() =
        a_scalar.gradient().getGrad() /
        (static_cast<ScalarType>(1) + a_scalar.value() * a_scalar.value());
    new_scalar.gradient().getHessian() =
        a_scalar.gradient().getHessian() /
            (static_cast<ScalarType>(1) + a_scalar.value() * a_scalar.value()) -
        a_scalar.gradient().getGrad() *
            (static_cast<ScalarType>(2) *
             a_scalar.gradient().getGrad().transpose() * a_scalar.value()) /
            ((static_cast<ScalarType>(1) +
              a_scalar.value() * a_scalar.value()) *
             (static_cast<ScalarType>(1) +
              a_scalar.value() * a_scalar.value()));
  } else {
    new_scalar.gradient() =
        a_scalar.gradient() /
        (static_cast<ScalarType>(1) + a_scalar.value() * a_scalar.value());
  }
  return new_scalar;
}

/*********************************** ATANH ***********************************/

template <class ScalarWithGradientType>
enable_if_t<!has_embedded_gradient<ScalarWithGradientType>::value,
            ScalarWithGradientType>
ArctanhMoments(const ScalarWithGradientType& a_scalar) {
  return atanh(a_scalar);
}

// template <>
// inline enable_if_t<!has_embedded_gradient<double>::value, double>
// ArctanhMoments(const double& a_scalar) {
//   return std::atanh(a_scalar);
// }
// template <>
// inline enable_if_t<!has_embedded_gradient<Quad_t>::value, Quad_t>
// ArctanhMoments(const Quad_t& a_scalar) {
//   return atanhq(a_scalar);
// }

template <class ScalarWithGradientType>
inline enable_if_t<has_embedded_gradient<ScalarWithGradientType>::value,
                   ScalarWithGradientType>
ArctanhMoments(const ScalarWithGradientType& a_scalar) {
  using ScalarType = typename ScalarWithGradientType::value_type;
  ScalarWithGradientType new_scalar(static_cast<ScalarType>(0));
  new_scalar.value() = atanh(a_scalar.value());
  if constexpr (ScalarWithGradientType::gradient_type::has_hessian) {
    new_scalar.gradient().getGrad() =
        a_scalar.gradient().getGrad() /
        (static_cast<ScalarType>(1) - a_scalar.value() * a_scalar.value());
    new_scalar.gradient().getHessian() =
        a_scalar.gradient().getHessian() /
            (static_cast<ScalarType>(1) - a_scalar.value() * a_scalar.value()) +
        a_scalar.gradient().getGrad() *
            (static_cast<ScalarType>(2) *
             a_scalar.gradient().getGrad().transpose() * a_scalar.value()) /
            ((static_cast<ScalarType>(1) -
              a_scalar.value() * a_scalar.value()) *
             (static_cast<ScalarType>(1) -
              a_scalar.value() * a_scalar.value()));
  } else {
    new_scalar.gradient() =
        a_scalar.gradient() /
        (static_cast<ScalarType>(1) - a_scalar.value() * a_scalar.value());
  }
  return new_scalar;
}

template <class GradientType>
using ScalarWithGradient = ScalarWithGradientBase<double, GradientType>;

}  // namespace IRL

#include "irl/geometry/general/scalar_with_gradient.tpp"

#endif