// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_GENERAL_MOMENTS_H_
#define IRL_MOMENTS_GENERAL_MOMENTS_H_

#include <ostream>

#include "irl/parameters/defined_types.h"

namespace IRL {

/// \brief A general moments class
/// that stores moments in row-major order.
/// For 3D, this is:
/// 1, x, y, z, x^2, xy, xz, y^2, yz, z^2, x^3, x^3 y, ...
/// For 2D, this is
/// 1, x, y, x^2, xy, y^2, x^3, x^3 y, ...
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM>
class GeneralMoments {
 public:
  static constexpr std::size_t linear_length =
      DIM == 3 ? (ORDER + 1) * (ORDER + 2) * (ORDER + 3) / 6
               : (ORDER + 1) * (ORDER + 2) / 2;
  using storage = std::array<double, linear_length>;

  /// \brief Default constructor.
  GeneralMoments(void);

  /// @brief  \breif Fill all moments with `a_value`
  static GeneralMoments fromScalarConstant(const double a_value);

  /// \brief Obtain GeneralMoments from the supplied geometry.
  template <class GeometryType>
  static GeneralMoments calculateMoments(GeometryType* a_geometry);

  /// \brief Return array of central moments.
  storage calculateCentralMoments(void) const;

  /// \brief Return scale and translation invariant moments,
  /// properly normalized by M^0
  storage calculateInvariantMoments(void) const;

  /// \brief Normalize currently stored moments by normalization used for
  /// scale invariance.
  void normalizeAsInvariant(void);

  /// \brief Return const reference to stored volume.
  double& volume(void);

  /// \brief Return const reference to stored volume.
  const double volume(void) const;

  /// \brief Multiply all moments by the volume (zeroeth moments)
  void multiplyByVolume(void);

  /// \brief Divide all moments by the volume (zeroeth moments)
  void normalizeByVolume(void);

  /// \brief Overload += operator to update moments.
  GeneralMoments& operator+=(const GeneralMoments& a_rhs);

  /// \brief Overload *= operator to multiply by constant double
  GeneralMoments& operator*=(const double a_rhs);

  /// \brief Overload /= operator to divide by constant double
  GeneralMoments& operator/=(const double a_rhs);

  /// \brief Overload assignment to assign constant value to moments.
  GeneralMoments& operator=(const double a_value);

  /// \brief Total number of moments.
  UnsignedIndex_t size(void) const;

  /// \brief Provide mutable access to underlying moments
  double& operator[](const UnsignedIndex_t index);

  /// \brief Provide const access to underlying moments
  double operator[](const UnsignedIndex_t index) const;

  /// \brief Default destructor.
  ~GeneralMoments(void) = default;

 private:
  storage moments_m;  // Stored moments in row-major order
};

template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM>
std::ostream& operator<<(std::ostream& out,
                         const GeneralMoments<ORDER, DIM>& a_volume);

/// \brief Overload + operator to add two geometric moments together
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM>
inline GeneralMoments<ORDER, DIM> operator+(
    const GeneralMoments<ORDER, DIM>& a_vm1,
    const GeneralMoments<ORDER, DIM>& a_vm2);
/// \brief Overload - operator to subtract one
/// geometric moment object from another.
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM>
inline GeneralMoments<ORDER, DIM> operator-(
    const GeneralMoments<ORDER, DIM>& a_vm1,
    const GeneralMoments<ORDER, DIM>& a_vm2);
/// \brief Overload * operator to multiply moments
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM>
inline GeneralMoments<ORDER, DIM> operator*(
    const double a_multiplier, const GeneralMoments<ORDER, DIM>& a_vm);
/// \brief Overload * operator to multiply moments
template <UnsignedIndex_t ORDER, UnsignedIndex_t DIM>
inline GeneralMoments<ORDER, DIM> operator*(
    const GeneralMoments<ORDER, DIM>& a_vm, const double a_multiplier);

template <UnsignedIndex_t ORDER>
using GeneralMoments3D = GeneralMoments<ORDER, 3>;

template <UnsignedIndex_t ORDER>
using GeneralMoments2D = GeneralMoments<ORDER, 2>;

}  // namespace IRL

#include "irl/moments/general_moments.tpp"

#endif  // IRL_MOMENTS_GENERAL_MOMENTS_H_ */
