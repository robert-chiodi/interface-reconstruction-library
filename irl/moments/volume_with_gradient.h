// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_MOMENTS_VOLUME_WITH_GRADIENT_H_
#define IRL_MOMENTS_VOLUME_WITH_GRADIENT_H_

#include <ostream>

#include "irl/moments/volume.h"
#include "irl/parameters/defined_types.h"

namespace IRL {

/// \brief A volume class which is just a double
/// with special properties. Allows more general
/// writing of functions.
template <class GradientType>
class VolumeWithGradient {
 public:
  using gradient_type = GradientType;

  /// \brief Default constructor.
  VolumeWithGradient(void);

  /// \brief Constructor that initializes volume to value,
  /// want to allow implicit casting of double to Volume here.
  constexpr VolumeWithGradient(const double a_value);

  constexpr VolumeWithGradient(const double a_volume,
                               const GradientType& a_gradient);

  static VolumeWithGradient fromScalarConstant(const double a_value);

  /// \brief Dummy function to allow general use along
  /// with VolumeMoments.
  void multiplyByVolume(void);

  /// \brief Dummy function to allow general use along
  /// with VolumeMoments.
  void normalizeByVolume(void);

  /// \brief Obtain un-normalized VolumeMoments from the supplied geometry.
  template <class GeometryType>
  static VolumeWithGradient calculateMoments(GeometryType* a_geometry);

  /// \brief Return value of stored volume.
  Volume& volume(void);

  /// \brief Return const reference to stored volume.
  const Volume& volume(void) const;

  /// \brief Return value of stored gradient.
  GradientType& gradient(void);

  /// \brief Return const reference to stored gradient.
  const GradientType& gradient(void) const;

  /// \brief Overload += operator to update volume.
  VolumeWithGradient& operator+=(const VolumeWithGradient& a_rhs);

  /// \brief Overload *= operator to multiply by constant double
  VolumeWithGradient& operator*=(const double a_rhs);

  /// \brief Overload /= operator to divide by constant double
  VolumeWithGradient& operator/=(const double a_rhs);

  /// \brief Allow implicit conversion to double.
  operator double() const;

  /// \brief Overload assignment to assign constant value to moments.
  VolumeWithGradient& operator=(const double a_value);

  /// \brief Default destructor.
  ~VolumeWithGradient(void) = default;

 private:
  Volume volume_m;          ///< \brief Volume of something
  GradientType gradient_m;  ///< \brief Gradient of the volume, with respect to
                            ///< some parameters
};

/// \brief Overload + operator to add two geometric moments together
template <class GradientType>
inline VolumeWithGradient<GradientType> operator+(
    const VolumeWithGradient<GradientType>& a_vm1,
    const VolumeWithGradient<GradientType>& a_vm2);
/// \brief Overload - operator to subtract one
/// geometric moment object from another.
template <class GradientType>
inline VolumeWithGradient<GradientType> operator-(
    const VolumeWithGradient<GradientType>& a_vm1,
    const VolumeWithGradient<GradientType>& a_vm2);
/// \brief Overload * operator to multiply volume/centroid
template <class GradientType>
inline VolumeWithGradient<GradientType> operator*(
    const double a_multiplier, const VolumeWithGradient<GradientType>& a_vm);
/// \brief Overload * operator to multiply volume/centroid
template <class GradientType>
inline VolumeWithGradient<GradientType> operator*(
    const VolumeWithGradient<GradientType>& a_vm, const double a_multiplier);

}  // namespace IRL

#include "irl/moments/volume_with_gradient.tpp"

#endif  // IRL_MOMENTS_VOLUME_WITH_GRADIENT_H_ */
