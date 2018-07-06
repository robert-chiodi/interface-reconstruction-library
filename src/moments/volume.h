// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_MOMENTS_VOLUME_H_
#define SRC_MOMENTS_VOLUME_H_

#include <ostream>

#include "src/parameters/defined_types.h"

namespace IRL {

/// \brief A volume class which is just a double
/// with special properties. Allows more general
/// writing of functions.
class Volume {
 public:
  /// \brief Default constructor.
  Volume(void);

  /// \brief Constructor that initializes volume to value,
  /// want to allow implicit casting of double to Volume here.
  constexpr Volume(const double a_value);

  static Volume fromScalarConstant(const double a_value);

  /// \brief Obtain Volume from the supplied geometry.
  template <class GeometryType>
  static Volume calculateMoments(GeometryType* a_geometry);

  /// \brief Dummy function to allow general use along
  /// with VolumeMoments.
  void multiplyByVolume(void);

  /// \brief Dummy function to allow general use along
  /// with VolumeMoments.
  void normalizeByVolume(void);

  /// \brief Overload += operator to update volume.
  Volume& operator+=(const Volume& a_rhs);

  /// \brief Overload *= operator to multiply by constant double
  Volume& operator*=(const double a_rhs);

  /// \brief Overload /= operator to divide by constant double
  Volume& operator/=(const double a_rhs);

  /// \brief Allow implicit conversion to double.
  operator double() const;

  /// \brief Overload assignment to assign constant value to moments.
  Volume& operator=(const double a_value);

  /// \brief Default destructor.
  ~Volume(void) = default;

 private:
  double volume_m;  ///< \brief Volume of something
};

std::ostream& operator<<(std::ostream& out, const Volume& a_volume);

}  // namespace IRL

#include "src/moments/volume.tpp"

#endif  // SRC_MOMENTS_VOLUME_H_ */
